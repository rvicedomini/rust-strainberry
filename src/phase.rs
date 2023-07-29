mod haplotree;
mod haplotype;
mod phasedblock;

use std::collections::VecDeque;
use std::fs;
use std::path::{Path,PathBuf};

use itertools::Itertools;
use rustc_hash::FxHashMap;

use rust_htslib::bam;
use rust_htslib::bam::Read;

use crate::cli::Options;
use crate::utils;
use crate::utils::BamRecordId;
use crate::utils::seq::{SeqInterval,SuccinctSeq};
use crate::variant::{Var,VarDict};

use phasedblock::PhasedBlock;

use self::haplotype::Haplotype;

// input: 
//   - read-alignment (bam)
//   - work-dir
//   - target sequences
//   - variant positions
//   - other options
// output:
//   - phased haplotypes (with read files)
//   - filtered variant positions

pub struct Phaser<'a> {
    bam: &'a Path,
    _target_sequences: &'a Vec<Vec<u8>>,
    work_dir: PathBuf,
    opts: &'a Options,
}

impl<'a> Phaser<'a> {

    pub fn new(bam: &'a Path, target_sequences: &'a Vec<Vec<u8>>, output_dir: &Path, opts: &'a Options) -> Phaser<'a> {

        let work_dir = output_dir.join("20-separate");
        if fs::create_dir_all(work_dir.as_path()).is_err() {
            eprintln!("Cannot create output directory: \"{}\"", work_dir.display());
            std::process::exit(1);
        };
        eprintln!("Phasing work directory: {}", work_dir.display());

        Phaser {
            bam,
            _target_sequences: target_sequences,
            work_dir,
            opts
        }
    }

    pub fn work_dir(&self) -> &Path {
        self.work_dir.as_path()
    }

    // TODO: split variants at positions that are too distant
    pub fn run(&self, variants: &VarDict) {

        for (&tid, target_variants) in variants {
            if target_variants.len() > 0 {
                let beg = target_variants.first().unwrap().pos;
                let end = target_variants.last().unwrap().pos + 1;
                let target_interval = SeqInterval{tid,beg,end};
                self.phase_interval(&target_interval, target_variants);
            }
        }

    }

    fn phase_interval(&self, target_interval:&SeqInterval, variants: &[Var]) {
        eprintln!("----- Phasing {target_interval}");
        let _haplotypes = self.phase_variants(target_interval,variants);
    }

    fn phase_variants(&self, target_interval:&SeqInterval, variants: &[Var]) -> Vec<Haplotype> {

        let mut haplotypes = vec![];

        if variants.len() == 0 {
            return haplotypes;
        }

        let mut succinct_records: FxHashMap<BamRecordId,SuccinctSeq> = FxHashMap::default();
        let mut lookback_positions: VecDeque<usize> = VecDeque::new();
        let mut var_iter = variants.iter();
        let mut phasedblock = PhasedBlock::new(target_interval.tid);

        let mut bam_reader = bam::IndexedReader::from_path(self.bam).unwrap();
        bam_reader.fetch((target_interval.tid as i32, target_interval.beg as i64, target_interval.end as i64)).expect(&format!("Cannot fetch target interval: {}",target_interval));

        let mut var = var_iter.next();

        for pileup in bam_reader.pileup().flatten() {
            
            if var.is_none() {
                break 
            };

            let target_pos = pileup.pos() as usize;
            
            let mut var_position = var.unwrap().pos;
            let mut var_nucleotides = var.unwrap().alleles.iter().map(|x| x.0 as u8).collect_vec();

            if var_position < target_pos {
                var = var_iter.next();
                if var.is_none() { break; }
                var_position = var.unwrap().pos;
                var_nucleotides = var.unwrap().alleles.iter().map(|x| x.0 as u8).collect_vec();
            }

            if target_pos < var_position {
                continue;
            }

            eprintln!("-----> SNV position = {var_position}");
            debug_assert_eq!(target_pos, var_position);

            // TODO: keep a list of "in-scope" succinct reads
            let edges = self.process_pileup(&pileup, &mut succinct_records, &var_nucleotides);
            lookback_positions.push_back(var_position);

            eprintln!("Nucleotides: {}, Edges: {:?}", String::from_utf8(var_nucleotides.clone()).unwrap(), edges.iter().map(|&(u,v)| (u as char, v as char)).collect_vec());

            // filter lookback positions too far away
            while var_position - lookback_positions[0] + 1 > self.opts.lookback {
                lookback_positions.pop_front();
            }

            debug_assert!(lookback_positions.len() > 0);

            // if this is the start a new phased block or there is no edge to extend current one
            if lookback_positions.len() == 1 || edges.len() == 0 {
                // FIXME? if there is no edge, why did I not save haplotypes ? Check old python code!!!
                phasedblock.init(var_position, var_nucleotides);
                continue
            }

            // eprintln!("Haplotypes:");
            // for ht in phasedblock.haplotypes().values() {
            //     eprintln!("{ht}");
            // }

            // identify haplotypes that cannot be extended with current edges
            let unsupported_haplotypes = phasedblock.haplotypes().iter()
                .filter(|&(_, ht)| edges.iter().all(|&(s,_)| s != ht.last_nuc() ))
                .map(|(&hid,_)| hid)
                .collect_vec();

            // possibly save haplotypes and start a new phased block
            if unsupported_haplotypes.len() > 0 && phasedblock.haplotypes().values().next().unwrap().size() >= 3 && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                eprintln!("Saving haplotypes after impossible extension");
                for ht in phasedblock.haplotypes().values() {
                    eprintln!("{ht}");
                }
                haplotypes.append(&mut phasedblock.drain());
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue
            }

            // remove haplotypes unsupported by the edges
            for hid in unsupported_haplotypes { 
                phasedblock.remove_haplotype(hid)
            }

            if phasedblock.haplotypes().len() < 2 {
                phasedblock.init(var_position, var_nucleotides);
                continue
            }

            let mut is_ambiguous = phasedblock.extend(var_position, edges);

            // discard haplotypes unsupported by the reads
            let back_pos = (var_position+1).checked_sub(self.opts.lookback).unwrap_or(0);
            let back_i = lookback_positions.partition_point(|&pos| pos < back_pos);
            let min_position = lookback_positions[back_i];
            
            // the following line could be slow -> consider the same approach in python version,
            // that is, use a "supporting_records" vector that stores the IDs of the last supporting succinct records
            // and to be updated within the process_pileup function
            succinct_records.retain(|_,sr_seq| sr_seq.positions().last().unwrap() == &var_position);
            
            let candidate_records = succinct_records.iter()
                .filter(|&(_,sr_seq)| sr_seq.positions()[0] <= min_position)
                .map(|(sr_id,_)| sr_id)
                .collect_vec();

            eprintln!("min_position: {}", min_position);
            eprintln!("supporting_records: {}", succinct_records.len());
            eprintln!("candidate_records: {}", candidate_records.len());

            let (unsupported_haplotypes, ambiguous_haplotypes) = self.validate_haplotypes(&succinct_records, &candidate_records, &phasedblock, min_position);

            if unsupported_haplotypes.len() > 0 && phasedblock.haplotypes().values().next().unwrap().size() >= 3 && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                phasedblock.split_and_init(0);
                // eprintln!("Saved haplotypes after finding unsupported ones:");
                // for ht in phasedblock.haplotypes().values() {
                //     eprintln!("{ht}");
                // }
                haplotypes.append(&mut phasedblock.drain());
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue;
            }

            for hid in unsupported_haplotypes {
                // eprintln!("unsupported: {}", phasedblock.haplotypes().get(&hid).unwrap());
                phasedblock.remove_haplotype(hid);
            }
            // for hid in ambiguous_haplotypes.iter() {
            //     eprintln!("ambiguous: {}", phasedblock.haplotypes().get(&hid).unwrap());
            // }
            is_ambiguous |= ambiguous_haplotypes.len() > 0;

            // if ambiguous extension, create a new phaseset (should a minimum of 3 SNVs be requested here too?)
            if is_ambiguous && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                let mut new_phasedblock = phasedblock.split_and_init(self.opts.lookback);
                haplotypes.append(&mut phasedblock.drain());
                std::mem::swap(&mut phasedblock, &mut new_phasedblock);
            } else { // discard ambiguous haplotypes if phased region was too short
                // eprintln!("Discarding short \"ambiguous\" haplotypes:");
                for hid in ambiguous_haplotypes {
                    // let ht = phasedblock.get(hid);
                    // eprintln!("{ht}");
                    phasedblock.remove_haplotype(hid);
                }
            }

            // eprintln!("Validated haplotypes:");
            // for ht in phasedblock.haplotypes().values() {
            //     eprintln!("{ht}");
            // }
        }

        // eprintln!("Saving remaining haplotypes:");
        // for ht in phasedblock.haplotypes().values() {
        //     eprintln!("{ht}");
        // }
        haplotypes.append(&mut phasedblock.drain());
        
        haplotypes
    }

    fn validate_haplotypes(&self, succinct_records: &FxHashMap<BamRecordId,SuccinctSeq>, candidate_records: &Vec<&BamRecordId>, phasedblock: &PhasedBlock, min_position: usize) -> (Vec<usize>,Vec<usize>) {

        fn get_best_haplotypes(sr: &SuccinctSeq, haplotypes:&Vec<&Haplotype>, min_position:usize) -> Vec<usize> {
            let mut sr_distances = vec![];
            for ht in haplotypes {
                let back_i = ht.raw_variants().partition_point(|snv| snv.pos < min_position);
                let back_pos = ht.raw_variants().get(back_i).unwrap().pos;
                let sr_left = sr.positions().partition_point(|pos| *pos < back_pos);
                let sr_right = sr.positions().partition_point(|pos| *pos <= ht.last_pos());
                let sr_nucleotides = sr.nucleotides().get(sr_left..sr_right).unwrap();
                let ht_nucleotides = ht.raw_variants().get(back_i..).unwrap();
                
                let hamming_dist = sr_nucleotides.iter().cloned()
                    .zip_eq(ht_nucleotides.iter().map(|snv| snv.nuc))
                    .filter(|(a,b)| a != b)
                    .count();

                if 3 * hamming_dist <= (sr_right-sr_left) {
                    sr_distances.push((ht.hid(), hamming_dist, back_i));
                }
            }
            if sr_distances.len() == 0 {
                return vec![]
            }
            let min_dist = sr_distances.iter().map(|(_,d,_)| *d).min().unwrap();
            sr_distances.into_iter()
                .filter(|(_,d,_)| *d == min_dist)
                .map(|(hid,_,_)| hid)
                .collect_vec()
        }
        
        let mut supporting: FxHashMap<usize,usize> = FxHashMap::default();
        let mut unambiguous: FxHashMap<usize,usize> = FxHashMap::default();
        let haplotypes = phasedblock.haplotypes().values().collect_vec();

        for &sr_id in candidate_records {
            let sr = &succinct_records[sr_id];
            let mut best_haplotypes = get_best_haplotypes(sr, &haplotypes, min_position);
            if best_haplotypes.len() == 1 {
                let hid: usize = best_haplotypes.pop().unwrap();
                unambiguous.entry(hid).and_modify(|cnt| *cnt+=1).or_insert(1);
                supporting.entry(hid).and_modify(|cnt| *cnt+=1).or_insert(1);
            }
            while let Some(hid) = best_haplotypes.pop() {
                supporting.entry(hid).and_modify(|cnt| *cnt+=1).or_insert(1);
            }
        }

        let unsupported_haplotypes = phasedblock.haplotypes().keys()
            .filter(|&hid| !supporting.contains_key(hid) || supporting[hid] < self.opts.min_alt_count)
            .cloned()
            .collect_vec();
        let ambiguous_haplotypes = phasedblock.haplotypes().keys()
            .filter(|&hid| supporting.contains_key(hid) && supporting[hid] >= self.opts.min_alt_count && (!unambiguous.contains_key(hid) || unambiguous[hid] < self.opts.min_alt_count))
            .cloned()
            .collect_vec();

        (unsupported_haplotypes, ambiguous_haplotypes)
    }


    fn process_pileup(&self, pileup: &bam::pileup::Pileup, succinct_records: &mut FxHashMap<BamRecordId,SuccinctSeq>, var_nucleotides: &Vec<u8>) -> Vec<(u8,u8)> {

        let position = pileup.pos() as usize;
        // let mut nuc_counter: FxHashMap<u8,usize> = FxHashMap::default();
        // let mut total_obs: usize = 0;

        let mut edge_total_obs: usize = 0;
        let mut edge_counter: FxHashMap<(u8,u8),usize> = FxHashMap::default();

        for alignment in pileup.alignments() {
            let record = alignment.record();

            if record.mapq() < self.opts.min_mapq
                || record.is_unmapped() 
                || record.is_secondary() 
                || record.is_quality_check_failed() 
                || record.is_duplicate() 
                || record.is_supplementary()
            {
                continue;
            }

            // total_obs += 1;
            
            let record_id = utils::bam_record_id(&record);
            let record_nuc = if !alignment.is_del() && !alignment.is_refskip() { 
                alignment.record().seq()[alignment.qpos().unwrap()] 
            } else { 
                b'-'
            };

            // *nuc_counter.entry(record_nuc).or_default() += 1;

            let srec = succinct_records.entry(record_id.clone())
                .or_insert(SuccinctSeq::build(&record_id.0, pileup.tid() as usize));
            srec.push(position, record_nuc);
            
            let srec_positions = srec.positions();
            let srec_len = srec.len();
            if srec_len > 1 && srec_positions[srec_len-1] - srec_positions[srec_len-2] + 1 <= self.opts.lookback {
                let srec_nucleotides = srec.nucleotides();
                let edge = (srec_nucleotides[srec_len-2], srec_nucleotides[srec_len-1]);
                if var_nucleotides.contains(&edge.1) {
                    edge_total_obs += 1;
                    edge_counter.entry(edge)
                        .and_modify(|cnt| *cnt += 1)
                        .or_insert(1);
                }
            }
        }

        // eprintln!("pileup at {}: {} >={}MAPQ", position, total_obs, self.opts.min_mapq);

        // let nucleotides = nuc_counter.iter()
        //     .filter(|&(_nuc,&cnt)| cnt >= self.opts.min_alt_count && (cnt as f64) >= self.opts.min_alt_frac * (total_obs as f64))
        //     .map(|(nuc,_cnt)| *nuc)
        //     .collect_vec();

        // let edge_total_obs = edge_counter.values().sum::<usize>() as f64;
        let x = edge_counter.iter().map(|(&(u,v),&cnt)| (format!("{}{}",u as char, v as char),cnt)).collect_vec();
        eprintln!("|E|:{edge_total_obs}, Edge: {x:?}");
        let edges = edge_counter.iter()
            .filter(|&(_edge,&cnt)| cnt >= self.opts.min_alt_count && (cnt as f64) >= self.opts.min_alt_frac * (edge_total_obs as f64))
            .map(|(&edge,_cnt)| edge)
            .collect_vec();

        edges
    }

}



// def _separate_reference_region(self, contig_id, contig_start, contig_end, variant_positions):

//     logger.debug(f'----- Separating {contig_id}:{contig_start}-{contig_end} -----')
//     contig_region_name = f'{contig_id}_{contig_start}-{contig_end}'

//     if len(variant_positions) == 0:
//         logger.debug(f'No SNV to be phased')
//         return {}, set()

//     logger.debug(f'Phasing haplotypes')
//     phaser = NextPhaser(self.bam, min_mapq=self.mapq, lookback=self.lookback, min_obs=self.min_obs, min_frac=self.min_frac)
//     haplotypes = phaser.phase_haplotypes(contig_id, variant_positions)
//     logger.debug(f'Phased haplotypes: {len(haplotypes)}')

//     logger.debug(f'Removing false haplotypes')
//     variant_positions = self.retrieve_variant_positions(haplotypes)
//     haplotypes = self.remove_false_haplotypes(haplotypes, variant_positions)
//     logger.debug(f'Remaining haplotypes: {len(haplotypes)}')

//     if len(haplotypes) == 0:
//         logger.debug(f'No haplotype remaining')
//         return {}, set()

//     logger.debug(f'Mapping reads to haplotypes')
//     haplotype_intervals = IntervalTree(Interval(ht.start(),ht.end(),ht) for ht in haplotypes)
//     segment_haplotypes, ambiguous_segments = local_separate_reads(phaser.succinct_reads.values(), haplotype_intervals, variant_positions)
//     logger.debug(f'Mapped={len(segment_haplotypes)} Ambiguous={len(ambiguous_segments)} Total={len(phaser.succinct_reads)}')

//     logger.debug(f'Haplotype extension')
//     hap_graph = HaploGraph(haplotypes, segment_haplotypes, self.lookback)
//     hap_graph.write_dot(os.path.join(self.outdir_dots,f'{contig_region_name}.haplotypes.partial.dot'))
//     haplotypes = hap_graph.scaffold_haplotypes(variant_positions, min_nreads=self.min_obs)
//     hap_graph.write_dot(os.path.join(self.outdir_dots,f'{contig_region_name}.haplotypes.final.dot'))
//     logger.debug(f'Extended haplotypes: {len(haplotypes)}')

//     logger.debug(f'Mapping reads to extended haplotypes')
//     haplotype_intervals = IntervalTree(Interval(ht.start(),ht.end(),ht) for ht in haplotypes)
//     segment_haplotypes, ambiguous_segments = local_separate_reads(phaser.succinct_reads.values(), haplotype_intervals, variant_positions)
//     logger.debug(f'Mapped={len(segment_haplotypes)} Ambiguous={len(ambiguous_segments)} Total={len(phaser.succinct_reads)}')

//     phasesets = defaultdict(list)
//     for ht in haplotypes:
//         phasesets[ht.region()].append(ht)

//     logger.debug(f'Writing reads')
//     for (ps_start,ps_end), ht_list in phasesets.items():
//         self.write_reads(contig_id, ps_start, ps_end, ht_list, segment_haplotypes)

//     return phasesets, variant_positions



// def separate_references(self):
//     for reference_id in self.reference_lengths:
//         reference_positions = self.variant_positions[reference_id]
//         reference_phasesets = dict()
//         reference_variant_positions = set()
//         for reference_start, reference_end in self.reference_intervals[reference_id]:
//             first = bisect.bisect_left(reference_positions, reference_start)
//             last = bisect.bisect_right(reference_positions, reference_end-1)
//             region_phasesets, region_positions = self._separate_reference_region(reference_id, reference_start, reference_end, reference_positions[first:last])
//             reference_phasesets |= region_phasesets
//             reference_variant_positions |= region_positions
//         self.phasesets[reference_id] = reference_phasesets
//         self.variant_positions[reference_id] = reference_variant_positions