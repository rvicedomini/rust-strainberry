use std::path::Path;
use std::sync::mpsc;
use std::thread;

use ahash::AHashMap as HashMap;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::htslib;

use crate::cli::Options;


// from ffforf: https://github.com/jguhlin/ffforf/blob/master/src/lib.rs
#[inline(always)]
fn complement(nuc: u8) -> u8 {
    if nuc != b'N' {
        if nuc & 2 != 0 {
            nuc ^ 4
        } else {
            nuc ^ 21
        }
    } else {
        nuc
    }
}

pub fn revcomp_inplace(seq: &mut [u8]) {
    seq.reverse();
    seq.iter_mut().for_each(|nuc| { *nuc = complement(*nuc) });
}


pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut rev_seq = seq.to_vec();
    revcomp_inplace(&mut rev_seq);
    rev_seq
}


pub fn load_bam_sequences(bam_path: &Path, read_index: &HashMap<String, usize>, opts: &Options) -> HashMap<usize,Vec<u8>> {
    
    let mut target_intervals = crate::bam::bam_target_intervals(bam_path);
    target_intervals.sort_unstable_by_key(|siv| siv.end - siv.beg);
    
    let (tx, rx) = mpsc::channel();
    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let target_intervals_ref = &target_intervals;
            scope.spawn(move || {
                let mut seq_dict = HashMap::new();
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for &siv in target_intervals_ref.iter().skip(thread_id).step_by(opts.nb_threads) {
                    bam_reader.fetch(siv.tid as u32).unwrap();
                    for record in bam_reader.rc_records().map(|x| x.expect("Failure parsing BAM file")) {
                        if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {
                            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                            let mut qseq = record.seq().as_bytes();
                            if record.is_reverse() {
                                revcomp_inplace(&mut qseq)
                            };
                            seq_dict.insert(read_index[&qname], qseq);
                        }
                    }
                }
                sender.send(seq_dict).unwrap();
            });
        }
    });

    (0..opts.nb_threads)
        .flat_map(|_| rx.recv().unwrap().into_iter())
        .collect()
}


pub fn build_read_index(bam_path: &Path) -> HashMap<String,usize> {

    let mut bam = IndexedReader::from_path(bam_path).unwrap();
    bam.fetch(".").unwrap();
    let read_index = bam.rc_records()
        .map(|x| x.expect("Error parsing BAM file"))
        .filter(|rec| rec.flags() & (htslib::BAM_FSECONDARY|htslib::BAM_FSUPPLEMENTARY) as u16 == 0) // consider only primary alignments
        .enumerate()
        .map(|(index,record)| (
            std::str::from_utf8(record.qname()).unwrap().to_string(),
            index
        )).collect();

    read_index
}
