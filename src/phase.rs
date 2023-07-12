mod haplotree;
mod haplotype;
mod phasedblock;

use std::collections::VecDeque;
use std::fs;
use std::path::{Path,PathBuf};

use crate::cli::Options;
use crate::utils::seq::SeqInterval;
use crate::variant::{Var,VarDict};

use phasedblock::PhasedBlock;

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
    target_sequences: &'a Vec<Vec<u8>>,
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

        Phaser { bam, target_sequences, work_dir, opts }
    }

    // TODO: split variants at positions that are too distant
    pub fn run(&mut self, variants: &VarDict) {

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
        eprintln!("------ Phasing {target_interval} -------");
        let haplotypes = self.phase_variants(target_interval,variants);
    }

    fn phase_variants(&self, target_interval:&SeqInterval, variants: &[Var]) {
        let mut lookback_positions: VecDeque<usize> = VecDeque::new();
        let mut var_iter = variants.iter();
        let mut var_position = var_iter.next().unwrap().pos;
        let mut phasedblock = PhasedBlock::new(target_interval.tid);
    }
    
}


// def _phase_snv_cluster(self, reference_id, reference_positions, supporting_sreads):
    
//     if len(reference_positions) == 0:
//         return self.haplotypes

//     logger.debug(f'## Phasing region: {reference_id}:{reference_positions[0]}-{reference_positions[-1]} ({len(reference_positions)} SNVs)')
    
//     lookback_positions = deque()
//     snv_position_iter = iter(reference_positions)
//     snv_position = next(snv_position_iter,None)
//     phased_block = PhasedBlock(reference_id)

//     with pysam.AlignmentFile(self.bamfile,'rb') as bam_handle:
//         for pileupcol in bam_handle.pileup(contig=reference_id, start=reference_positions[0], stop=reference_positions[-1]+1, stepper='all', min_base_quality=0):
            
//             # possibly update snv_position
//             if snv_position is None:
//                 break
//             if snv_position < pileupcol.reference_pos:
//                 snv_position = next(snv_position_iter,None)
//                 if snv_position is None:
//                     break
//             if pileupcol.reference_pos < snv_position:
//                 continue
            
//             # logger.debug(f'-----> SNV position = {snv_position}')
//             assert(pileupcol.reference_pos == snv_position)

//             # load pileup column
//             nucleotides, edges, supporting_sreads = self._process_column_pileup(pileupcol,supporting_sreads)
//             lookback_positions.append(snv_position)

//             # while snv_position-lookback_positions[0]+1 > self.lookback:
//             #     lookback_positions.popleft()
//             # assert(len(lookback_positions) > 0)

//             # start of a new phaseset or no available edge to extend current one
//             if len(lookback_positions) == 1 or len(edges) == 0:
//                 phased_block.init_haplotypes(snv_position, nucleotides)
//                 continue

//             # remove/save haplotypes that cannot be extended
//             unsupported_haplotypes = [ht_id for ht_id in phased_block.haplotypes if phased_block.haplotypes[ht_id].sequence[-1] not in map(operator.itemgetter(0),edges)]
//             if len(unsupported_haplotypes) > 0 and phased_block.haplotypes[unsupported_haplotypes[0]].size() >= 3 and (snv_position-phased_block.hap_begin+1 > self.lookback):
//                 # logger.debug(f'Saved haplotypes after impossible extension:')
//                 # for ht_id, ht in phased_block.haplotypes.items():
//                 #     logger.debug(f'{ht_id}\t{ht.seqstr()}\t[{ht.positions[ht.start_idx]}-{ht.positions[-1]}]')
//                 self._save_haplotypes(phased_block.haplotypes.values())
//                 phased_block.init_haplotypes(snv_position, nucleotides)
//                 lookback_positions = deque([snv_position])
//                 continue
//             for ht_id in unsupported_haplotypes:
//                 phased_block.remove_haplotype(ht_id)

//             if len(phased_block.haplotypes) < 2:
//                 phased_block.init_haplotypes(snv_position, nucleotides)
//                 continue

//             # extend viable haplotypes
//             is_ambiguous = phased_block.extend(snv_position, edges)

//             # discard unsupported haplotypes
//             back_i = bisect.bisect_left(lookback_positions,max(0,snv_position-self.lookback+1))
//             min_position = lookback_positions[back_i]
//             supporting_sreads = [sread for sread in supporting_sreads if sread.positions[-1] == snv_position]
//             # logger.debug(f'supporting_sreads = {len(supporting_sreads)}')
//             # logger.debug(f'min_position = {min_position}')
//             # logger.debug(f'candidate_reads = {len(list((sread for sread in supporting_sreads if sread.positions[0] <= min_position)))}')
//             candidate_reads = [sread for sread in supporting_sreads if sread.positions[0] <= min_position]
//             unsupported_haplotypes, ambiguous_haplotypes = self._validate_haplotypes(candidate_reads, phased_block.haplotypes, min_position)
//             # logger.debug(f'  haplotypes: candidates={len(phased_block.haplotypes)} unsupported={len(unsupported_haplotypes)} ambiguous={len(ambiguous_haplotypes)}')
            
//             if len(unsupported_haplotypes) > 0 and phased_block.haplotypes[unsupported_haplotypes[0]].size() >= 3 and (snv_position-phased_block.hap_begin+1 > self.lookback):
//                 phased_block.split_and_init(snv_position)
//                 # logger.debug(f'Saved haplotypes after finding unsupported ones:')
//                 # for ht_id, ht in phased_block.haplotypes.items():
//                 #     logger.debug(f'{ht_id}\t{ht.seqstr()}\t[{ht.positions[ht.start_idx]}-{ht.positions[-1]}]')
//                 self._save_haplotypes(phased_block.haplotypes.values())
//                 phased_block.init_haplotypes(snv_position, nucleotides)
//                 lookback_positions = deque([snv_position])
//                 continue

//             for ht_id in unsupported_haplotypes:
//                 # logger.debug(f'unsupported: {phased_block.haplotypes[ht_id].seqstr()} [{phased_block.haplotypes[ht_id].positions[0]}-{phased_block.haplotypes[ht_id].positions[-1]}]')
//                 phased_block.remove_haplotype(ht_id)
//             # for ht_id in ambiguous_haplotypes:
//             #     logger.debug(f'ambiguous: [{phased_block.haplotypes[ht_id].positions[0]}-{phased_block.haplotypes[ht_id].positions[-1]}] {phased_block.haplotypes[ht_id].sequence}')
//             is_ambiguous = is_ambiguous or (len(ambiguous_haplotypes) > 0)

//             # if ambiguous extension, create new phaseset (should a minimum number of 3 SNVs be requested here too?)
//             if is_ambiguous and (snv_position-phased_block.hap_begin+1 > self.lookback):
//                 # logger.debug(f'ambiguous extension, create new phased_block (snv_pos={snv_position} hap_beg={phased_block.hap_begin} lb={self.lookback})')
//                 # logger.debug(f'ambiguous haplotypes={ambiguous_haplotypes} edges={edges}')
//                 new_phased_block = phased_block.split_and_init(snv_position, self.lookback)
//                 # logger.debug(f'Saved haplotypes at ambiguous extension:')
//                 # for ht_id, ht in phased_block.haplotypes.items():
//                 #     logger.debug(f'{ht_id}\t{ht.seqstr()}\t[{ht.positions[ht.start_idx]}-{ht.positions[-1]}]')
//                 self._save_haplotypes(phased_block.haplotypes.values())
//                 phased_block = new_phased_block
//             elif len(ambiguous_haplotypes) > 0: # discard ambiguous haplotypes if phased region was too short
//                 # logger.debug(f'Discarding short "ambiguous" haplotypes:')
//                 for ht_id in ambiguous_haplotypes:
//                     ht = phased_block.haplotypes[ht_id]
//                     # logger.debug(f'{ht_id}\t{ht.seqstr()}\t{ht.positions}')
//                     phased_block.remove_haplotype(ht_id)

//             # logger.debug(f'Validated haplotypes:')
//             # for ht_id, ht in phased_block.haplotypes.items():
//             #     logger.debug(f'{ht_id}\t{ht.seqstr()}\t{ht.positions}')
        
//     # save remaning haplotypes
//     # logger.debug(f'Saving remaining haplotypes:')
//     # for ht_id, ht in phased_block.haplotypes.items():
//     #     logger.debug(f'{ht_id}\t{ht.seqstr()}\t[{ht.positions[ht.start_idx]}-{ht.positions[-1]}]')
//     self._save_haplotypes(phased_block.haplotypes.values())
//     return supporting_sreads


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