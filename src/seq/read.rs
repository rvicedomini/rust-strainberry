use std::path::Path;
use std::sync::mpsc;
use std::thread;

use itertools::Itertools;
use rust_htslib::bam::{Read,IndexedReader};
use rustc_hash::FxHashMap;

use crate::cli::Options;
use crate::utils;


#[inline(always)]
fn complement(nuc: u8) -> u8 {
    match nuc {
        b'a'|b'A' => b'T',
        b'c'|b'C' => b'G',
        b't'|b'T' => b'A',
        b'g'|b'G' => b'C',
        _ => b'N'
    }
}

#[inline(always)]
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rev_seq = Vec::with_capacity(seq.len());
    rev_seq.extend(seq.iter().rev().map(|&nuc| complement(nuc)));
    rev_seq
}


pub fn load_bam_sequences(bam_path: &Path, opts: &Options) -> FxHashMap<String,Vec<u8>> {
    
    let target_intervals = utils::bam_target_intervals(bam_path);
    let target_intervals = target_intervals.into_iter()
        .sorted_unstable_by_key(|siv| siv.end - siv.beg)
        .collect_vec();
    
    let (tx, rx) = mpsc::channel();
    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let target_intervals_ref = &target_intervals;
            scope.spawn(move || {
                let mut seq_dict = FxHashMap::default();
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for &siv in target_intervals_ref.iter().skip(thread_id).step_by(opts.nb_threads) {
                    bam_reader.fetch(siv.tid as u32).unwrap();
                    for record in bam_reader.rc_records().map(|x| x.expect("Failure parsing BAM file")) {
                        if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {
                            let qname = String::from_utf8_lossy(record.qname()).to_string();
                            let mut qseq = record.seq().as_bytes();
                            if record.is_reverse() {
                                qseq = reverse_complement(&qseq)
                            };
                            seq_dict.insert(qname, qseq);
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

