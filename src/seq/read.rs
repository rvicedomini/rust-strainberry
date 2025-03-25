use std::path::Path;
use std::sync::mpsc;
use std::thread;

use ahash::AHashMap as HashMap;
use rust_htslib::bam::{Read,IndexedReader};

use crate::cli::Options;
use crate::utils;


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

fn revcomp_inplace(seq: &mut [u8]) {
    seq.reverse();
    seq.iter_mut().for_each(|nuc| { *nuc = complement(*nuc) });
}


fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut rev_seq = seq.to_vec();
    revcomp_inplace(&mut rev_seq);
    rev_seq
}


pub fn load_bam_sequences(bam_path: &Path, opts: &Options) -> HashMap<String,Vec<u8>> {
    
    let mut target_intervals = utils::bam_target_intervals(bam_path);
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
                            let qname = String::from_utf8_lossy(record.qname()).to_string();
                            let mut qseq = record.seq().as_bytes();
                            if record.is_reverse() {
                                revcomp_inplace(&mut qseq)
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

