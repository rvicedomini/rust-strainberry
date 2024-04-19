pub mod haplograph;
pub mod haplotree;
pub mod haplotype;
pub mod phasedblock;

use std::collections::VecDeque;
use std::fs;
use std::sync::mpsc;
use std::thread;
use std::path::{Path,PathBuf};

use itertools::Itertools;
use rustc_hash::{FxHashSet,FxHashMap};

use rust_htslib::bam;
use rust_htslib::bam::Read;

use crate::cli::Options;
use crate::utils;
use crate::utils::BamRecordId;
use crate::seq::{SeqInterval,SuccinctSeq};
use crate::variant::{Var,VarDict};

use phasedblock::PhasedBlock;

// use self::haplograph::HaploGraph;
use self::haplotype::{Haplotype, HaplotypeId};

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
    bam_path: &'a Path,
    target_names: &'a Vec<String>,
    target_intervals: &'a Vec<SeqInterval>,
    read_sequences: &'a FxHashMap<String,Vec<u8>>,
    work_dir: PathBuf,
    opts: &'a Options,
}

impl<'a> Phaser<'a> {

    pub fn new(bam_path: &'a Path, target_names: &'a Vec<String>, target_intervals: &'a Vec<SeqInterval>, read_sequences: &'a FxHashMap<String, Vec<u8>>, output_dir: &Path, opts: &'a Options) -> Phaser<'a> {

        let work_dir = output_dir.join("20-separate");
        if fs::create_dir_all(work_dir.as_path()).is_err() {
            eprintln!("Cannot create output directory: \"{}\"", work_dir.display());
            std::process::exit(1);
        };
        eprintln!("Phasing work directory: {}", work_dir.display());

        Phaser {
            bam_path,
            target_names,
            target_intervals,
            read_sequences,
            work_dir,
            opts
        }
    }


    pub fn work_dir(&self) -> &Path {
        self.work_dir.as_path()
    }

    
    pub fn phase(&self, variants: &VarDict) -> FxHashMap<HaplotypeId,Haplotype> {

        let mut target_variants = vec![];
        for siv in self.target_intervals {
            if let Some(vars) = variants.get(&siv.tid) {
                let left = vars.partition_point(|v| v.pos < siv.beg);
                let right = vars.partition_point(|v: &Var| v.pos < siv.end);
                if right - left >= self.opts.min_snv {
                    target_variants.push((siv, &vars[left..right]));
                }
            }
        }
        target_variants.sort_unstable_by_key(|(_,vars)| vars.len());

        let nb_threads = std::cmp::min(self.opts.nb_threads, target_variants.len());
        let (tx, rx) = mpsc::channel();
        thread::scope(|scope| {
            for thread_id in 0..nb_threads {
                let sender = tx.clone();
                let target_variants_ref = &target_variants;
                scope.spawn(move || {
                    for &(target_interval, variants) in target_variants_ref[thread_id..].iter().step_by(nb_threads) {
                        let haplotypes = self.phase_interval(target_interval, variants);
                        sender.send(haplotypes).unwrap();
                    }
                });
            }
        });

        (0..target_variants.len())
            .flat_map(|_| rx.recv().unwrap())
            .collect()
    }

    fn phase_interval(&self, target_interval:&SeqInterval, variants: &[Var]) -> FxHashMap<HaplotypeId,Haplotype> {
        // eprintln!("  -> Phasing {target_interval} ({} variants)", variants.len());
        let (haplotypes, succinct_records) = self.phase_variants(target_interval,variants);
        // eprintln!("  -> Haplotypes in {target_interval}: {}", haplotypes.len());
        let variant_positions = self.variant_positions(&haplotypes);
        let haplotypes = self.remove_false_haplotypes(haplotypes, &variant_positions);
        // eprintln!("  -> Filtered Haplotypes in {target_interval}: {}", haplotypes.len());
        
        // let haplotype_intervals = IntervalTree::from_iter(
        //     haplotypes.values().map(|ht| (ht.beg()..ht.end(), ht))
        // );

        let (sread_haplotypes, _ambiguous_reads) = self::separate_reads(&succinct_records, &haplotypes, &variant_positions, self.opts.min_shared_snv);

        // TODO: Merge contiguous haplotypes when it is not ambiguous to do so
        // let mut haplograph = HaploGraph::new(&haplotypes, &sread_haplotypes);
        // let dot_file = format!("{}_{}-{}.dot", target_interval.tid, target_interval.beg, target_interval.end);
        // let dot_path = self.work_dir().join(dot_file.as_str());
        // let _ = haplograph.write_dot(&dot_path);

        // let mut phasesets: FxHashMap<SeqInterval,Vec<Haplotype>> = FxHashMap::default();
        // for ht in haplotypes.into_values() {
        //     phasesets.entry(ht.region())
        //         .or_default()
        //         .push(ht);
        // }
        
        self.write_reads(&haplotypes, &sread_haplotypes).unwrap();
        
        haplotypes
    }


    fn write_reads(&self, haplotypes: &FxHashMap<HaplotypeId,Haplotype>, sread_haplotypes: &FxHashMap<BamRecordId,Vec<HaplotypeId>>) -> std::io::Result<()> {
        
        let mut read_files = FxHashMap::default();
        
        for ht in haplotypes.values() {
            let ht_id = ht.uid();
            let ht_file_gz = format!("{}_{}-{}_h{}.fa.gz", self.target_names[ht_id.tid], ht_id.beg, ht_id.end, ht_id.hid);
            let ht_file_path = self.work_dir().join(ht_file_gz.as_str());
            read_files.insert(ht_id, utils::get_file_writer(&ht_file_path));
        }
        
        for (record_id, ht_list) in sread_haplotypes.iter() {
            if let Some(ht_id) = ht_list.iter().find(|&ht_id| read_files.contains_key(ht_id)) {
                let out = read_files.get_mut(ht_id).unwrap();
                out.write_all(format!(">{}\n", &record_id.0).as_bytes())?;
                out.write_all(self.read_sequences[&record_id.0].as_slice())?;
                out.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    fn remove_false_haplotypes(&self, mut haplotypes: FxHashMap<HaplotypeId,Haplotype>, positions: &FxHashSet<(usize,usize)>) -> FxHashMap<HaplotypeId,Haplotype> {
        haplotypes.retain(|_, ht| {
            ht.variants()
                .iter()
                .map(|snv| (ht.tid(),snv.pos))
                .any(|pos| positions.contains(&pos))
        });
        haplotypes
    }

    fn variant_positions(&self, haplotypes: &FxHashMap<HaplotypeId,Haplotype>) -> FxHashSet<(usize,usize)> {
        
        // let haplo_nodes: FxHashSet<&Snv> = haplotypes.values()
        //     .flat_map(|ht| ht.variants())
        //     .collect();

        let mut counter: FxHashMap<_,usize> = FxHashMap::default();
        for (tid,snv_pos) in haplotypes.values().flat_map(|ht| ht.variants().iter().map(|snv| (ht.tid(),snv.pos))) {
            counter.entry((tid,snv_pos))
                .and_modify(|cnt| *cnt+=1)
                .or_insert(1);
        }

        counter.into_iter()
            .filter(|&(_,cnt)| cnt > 1)
            .map(|(pos,_)| pos)
            .collect()
    }

    pub fn phase_variants(&self, target_interval:&SeqInterval, variants: &[Var]) -> (FxHashMap<HaplotypeId,Haplotype>,Vec<SuccinctSeq>)  {

        if variants.is_empty() {
            return (FxHashMap::default(),vec![]);
        }

        let mut haplotypes = FxHashMap::default();
        let mut succinct_records: FxHashMap<BamRecordId,SuccinctSeq> = FxHashMap::default();

        let mut supporting_sreads: FxHashSet<BamRecordId> = FxHashSet::default();
        let mut lookback_positions: VecDeque<usize> = VecDeque::new();
        let mut var_iter = variants.iter();
        let mut phasedblock = PhasedBlock::new(target_interval.tid);

        let mut bam_reader = bam::IndexedReader::from_path(self.bam_path).unwrap();
        bam_reader.fetch((target_interval.tid as i32, target_interval.beg as i64, target_interval.end as i64)).unwrap_or_else(|_| panic!("Cannot fetch target interval: {}",target_interval));

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

            // eprintln!("-----> SNV position = {var_position}");
            debug_assert_eq!(target_pos, var_position);

            // TODO: keep a list of "in-scope" succinct reads
            let edges = self.process_pileup(&pileup, &mut succinct_records, &mut supporting_sreads, &var_nucleotides);
            lookback_positions.push_back(var_position);

            // eprintln!("Nucleotides: {}, Edges: {:?}", String::from_utf8(var_nucleotides.clone()).unwrap(), edges.iter().map(|&(u,v)| (u as char, v as char)).collect_vec());

            // filter lookback positions too far away
            // while var_position - lookback_positions[0] + 1 > self.opts.lookback {
            //     lookback_positions.pop_front();
            // }

            debug_assert!(!lookback_positions.is_empty());

            // if this is the start a new phased block
            if lookback_positions.len() == 1 {
                phasedblock.init(var_position, var_nucleotides);
                continue
            }

            // if no edge within lookback bases, possibly save haplotypes
            if edges.is_empty() && phasedblock.haplotypes().values().next().unwrap().raw_size() >= self.opts.min_snv && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                for ht in phasedblock.drain() {
                    haplotypes.insert(ht.uid(), ht);
                }
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
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
            if !unsupported_haplotypes.is_empty() && phasedblock.haplotypes().values().next().unwrap().raw_size() >= self.opts.min_snv && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                // eprintln!("Saving haplotypes after impossible extension");
                // for ht in phasedblock.haplotypes().values() {
                //     eprintln!("{ht}");
                // }
                for ht in phasedblock.drain() {
                    haplotypes.insert(ht.uid(), ht);
                }
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
            let back_pos = (var_position+1).saturating_sub(self.opts.lookback);
            let back_i = lookback_positions.partition_point(|&pos| pos < back_pos);
            let min_position = lookback_positions[back_i];
            
            supporting_sreads.retain(|sr_id| succinct_records[sr_id].positions().last().unwrap() == &var_position);
            let candidate_records = supporting_sreads.iter()
                .filter(|&sr_id| succinct_records[sr_id].positions()[0] <= min_position)
                .collect_vec();
            
            // let candidate_records = succinct_records.iter()
            //     .filter(|&(_,sr_seq)| sr_seq.positions()[0] <= min_position)
            //     .map(|(sr_id,_)| sr_id)
            //     .collect_vec();

            // eprintln!("min_position: {}", min_position);
            // eprintln!("supporting_records: {}", succinct_records.len());
            // eprintln!("candidate_records: {}", candidate_records.len());

            let (unsupported_haplotypes, ambiguous_haplotypes) = self.validate_haplotypes(&succinct_records, &candidate_records, &phasedblock, min_position);

            if !unsupported_haplotypes.is_empty() && phasedblock.haplotypes().values().next().unwrap().raw_size() >= self.opts.min_snv && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                phasedblock.split_and_init(0);
                // eprintln!("Saved haplotypes after finding unsupported ones:");
                // for ht in phasedblock.haplotypes().values() {
                //     eprintln!("{ht}");
                // }
                for ht in phasedblock.drain() {
                    haplotypes.insert(ht.uid(), ht);
                }
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
            is_ambiguous |= !ambiguous_haplotypes.is_empty();

            // if ambiguous extension, create a new phaseset (should a minimum of 3 SNVs be requested here too?)
            if is_ambiguous && phasedblock.haplotypes().values().next().unwrap().raw_size() >= self.opts.min_snv && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                let mut new_phasedblock = phasedblock.split_and_init(self.opts.lookback);
                for ht in phasedblock.drain() {
                    haplotypes.insert(ht.uid(), ht);
                }
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
        if phasedblock.haplotypes().values().next().unwrap().raw_size() >= self.opts.min_snv {
            for ht in phasedblock.drain() {
                haplotypes.insert(ht.uid(), ht);
            }
        }
        
        (haplotypes, succinct_records.into_values().collect_vec())
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
            if sr_distances.is_empty() {
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

    fn process_pileup(&self, pileup: &bam::pileup::Pileup, succinct_records: &mut FxHashMap<BamRecordId,SuccinctSeq>, supporting_sreads: &mut FxHashSet<BamRecordId>, var_nucleotides: &[u8]) -> Vec<(u8,u8)> {

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

            if !succinct_records.contains_key(&record_id) {
                let srec = SuccinctSeq::build(record_id.clone(), pileup.tid() as usize);
                succinct_records.insert(record_id.clone(), srec);
                supporting_sreads.insert(record_id.clone());
            }

            let srec = succinct_records.get_mut(&record_id).unwrap();
            srec.push(position, record_nuc);

            let srec_positions = srec.positions();
            let srec_len = srec.len();
            if srec_len > 1 && srec_positions[srec_len-1] - srec_positions[srec_len-2] < self.opts.lookback {
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
        let _x = edge_counter.iter().map(|(&(u,v),&cnt)| (format!("{}{}",u as char, v as char),cnt)).collect_vec();
        // eprintln!("|E|:{edge_total_obs}, Edge: {_x:?}");
        let edges = edge_counter.iter()
            .filter(|&(_edge,&cnt)| cnt >= self.opts.min_alt_count && (cnt as f64) >= self.opts.min_alt_frac * (edge_total_obs as f64))
            .map(|(&edge,_cnt)| edge)
            .collect_vec();

        edges
    }

}


pub fn separate_reads(succinct_records: &[SuccinctSeq], haplotypes: &FxHashMap<HaplotypeId,Haplotype>, variant_positions: &FxHashSet<(usize,usize)>, min_shared_snv: usize)
    ->  (FxHashMap<BamRecordId,Vec<HaplotypeId>>, FxHashSet<String>) {

    let haplotypes = haplotypes.values().sorted_unstable_by_key(|ht| ht.uid()).collect_vec();
    let mut sread_haplotypes = FxHashMap::default();
    let mut ambiguous = FxHashSet::default();
    for sread in succinct_records {
        let best_hits = best_sread_haplotypes(sread, &haplotypes, variant_positions, min_shared_snv);
        if best_hits.iter().any(|(_hid,_ht_dist,alt_hits)| *alt_hits > 1) {
            ambiguous.insert(sread.record_id().0);
            continue
        }
        sread_haplotypes.entry(sread.record_id())
            .or_insert(vec![])
            .extend(best_hits.into_iter().map(|(hid,_,_)| hid));
    }
    (sread_haplotypes,ambiguous)
}

fn best_sread_haplotypes(sread: &SuccinctSeq, haplotypes: &[&Haplotype], variant_positions: &FxHashSet<(usize,usize)>, min_shared_snv: usize)
    -> Vec<(HaplotypeId,usize,usize)> {
    
    let mut candidates = FxHashMap::default();
    let sread_positions: FxHashSet<(usize,usize)> = sread.positions().iter().map(|var_pos| (sread.tid(),*var_pos)).collect();
    
    let idx = haplotypes.partition_point(|ht| ht.tid() < sread.tid() || (ht.tid() == sread.tid() && ht.end() <= sread.beg()));
    for ht in haplotypes[idx..].iter().take_while(|ht| ht.beg() < sread.end()) {
        let ht_left = ht.variants().partition_point(|snv| snv.pos < *sread.positions().first().unwrap());
        let ht_right = ht.variants().partition_point(|snv| snv.pos < *sread.positions().last().unwrap());
        let ht_positions: FxHashSet<(usize,usize)> = ht.variants()[ht_left..ht_right].iter().map(|snv| (ht.tid(),snv.pos)).collect();
        let shared_positions = ht_positions.intersection(&sread_positions).cloned().collect();
        if variant_positions.intersection(&shared_positions).count() >= min_shared_snv {
            let sread_nucleotides = sread.positions().iter()
                .enumerate()
                .filter(|&(_i,pos)| shared_positions.contains(&(sread.tid(),*pos)))
                .map(|(i,_)| *sread.nucleotides().get(i).unwrap())
                .collect_vec();
            let ht_nucleotides = ht.variants().iter()
                .filter(|&snv| shared_positions.contains(&(ht.tid(),snv.pos)))
                .map(|snv| snv.nuc)
                .collect_vec();
            let hamming_dist = sread_nucleotides.into_iter()
                .zip_eq(ht_nucleotides.iter())
                .filter(|&(a,b)| a != *b)
                .count();
            // hamming distance should not exceed 1/3 of the sequence length
            if 3 * hamming_dist <= ht_nucleotides.len() {
                let range = ht.beg()..ht.end();
                candidates.entry(range)
                    .or_insert(vec![])
                    .push((ht.uid(),hamming_dist));
            }
        }
    }

    let mut best_hits = vec![];
    for hits in candidates.into_values() {
        let best_dist = hits.iter().min_by_key(|&(_,ht_dist)| ht_dist).unwrap().1;
        let region_best_hits = hits.iter()
            .filter(|&(_,ht_dist)| *ht_dist == best_dist)
            .cloned()
            .collect_vec();
        // println!("sread:{}, dist:{}, hits:{:?}", sread.name(), best_dist, region_best_hits);
        best_hits.push((
            region_best_hits.last().unwrap().0,
            region_best_hits.last().unwrap().1,
            region_best_hits.len()
        ));
    }
    best_hits
}

