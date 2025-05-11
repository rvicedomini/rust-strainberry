mod alignment;
mod polish;
mod sequence;
mod window;

use ahash::AHashMap as HashMap;
use itertools::Itertools;
use rust_htslib::bam::record::CigarString;

use crate::awarecontig::{AwareAlignment, AwareContig};
use crate::cli::Options;
use crate::phase::haplotype::{HaplotypeId, Haplotype};
use crate::seq::SeqDatabase;

pub fn build_haplotigs(ref_db: &SeqDatabase, read_db: &SeqDatabase, haplotypes: &HashMap<HaplotypeId,Haplotype>, aware_contigs: &[AwareContig], read2aware: &HashMap<usize,Vec<AwareAlignment>>, opts: &Options) -> HashMap<HaplotypeId,Vec<u8>> {

    rayon::ThreadPoolBuilder::new().num_threads(opts.nb_threads).build_global().unwrap();

    let haplotypes = haplotypes.keys().cloned().sorted().collect_vec();
    let mut hap_index = HashMap::new();
    let hap_sequences = haplotypes.iter().enumerate().map(|(idx, hid)| {
        hap_index.insert(*hid, idx);
        ref_db.sequences[hid.tid][hid.beg..hid.end].to_vec()
    }).collect_vec();

    let mut alignments = Vec::new();
    for read_alignments in read2aware.values() {
        for a in read_alignments {
            let ctg = &aware_contigs[a.aware_id];
            if a.is_ambiguous || !ctg.is_phased() {
                continue
            }

            let hid = aware_contigs[a.aware_id].haplotype_id().unwrap();
            let hid_idx = hap_index[&hid];
            let hid_sequence = &hap_sequences[hid_idx];

            assert!(a.target_beg <= a.target_end);
            assert!(hid.beg <= a.target_beg && a.target_end <= hid.end);

            let length = std::cmp::max(a.query_end - a.query_beg, a.target_end-a.target_beg);

            alignments.push( alignment::Alignment {
                query_idx: a.query_idx,
                query_len: a.query_len,
                query_beg: a.query_beg,
                query_end: a.query_end,
                strand: a.strand,
                target_idx: hid_idx,
                target_len: hid_sequence.len(),
                target_beg: a.target_beg - ctg.beg(),
                target_end: a.target_end - ctg.beg(),
                matches: 0,         // dummy, not supposed to be used
                mapping_length: 0,  // dummy, not supposed to be used
                mapq: 60,           // dummy, not supposed to be used
                length,
                identity: 0.0,      // dummy, not supposed to be used
                cigar: CigarString(Vec::new()),
                breaking_points: Vec::new()
            });
        }
    }

    let mut polisher = polish::Polisher::new(&hap_sequences, &read_db.sequences, alignments);

    spdlog::debug!("initializing polisher");
    polisher.initialize();

    spdlog::debug!("polishing windows");
    let mut polished_sequences = polisher.polish();

    hap_index.drain()
        .map(|(hid, hap_idx)| (hid, std::mem::take(&mut polished_sequences[hap_idx])))
        .collect()
}
