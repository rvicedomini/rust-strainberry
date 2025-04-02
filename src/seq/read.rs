use std::path::Path;
use std::thread;

use ahash::AHashMap as HashMap;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::htslib;

use crate::cli::Options;
use crate::seq::bitseq::BitSeq;

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


pub fn load_bam_sequences(bam_path: &Path, read_index: &HashMap<String, usize>, opts: &Options) -> Vec<BitSeq> {
    
    let target_intervals = crate::bam::bam_target_intervals(bam_path);
    // target_intervals.sort_unstable_by_key(|siv| siv.end - siv.beg);

    let mut read_sequences = Vec::new();
    read_sequences.resize_with(read_index.len(), BitSeq::default);
    let read_sequences_slice = crate::utils::UnsafeSlice::new(&mut read_sequences);

    thread::scope(|scope| {
        let mut thread_handles = Vec::with_capacity(opts.nb_threads);
        for thread_id in 0..opts.nb_threads {
            let target_intervals_ref = &target_intervals;
            let read_sequences_ref = &read_sequences_slice;
            let handle = scope.spawn(move || {
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for &siv in target_intervals_ref.iter().skip(thread_id).step_by(opts.nb_threads) {
                    bam_reader.fetch(siv.tid as u32).unwrap();
                    for record in bam_reader.rc_records().map(|x| x.expect("Failure parsing BAM file")) {
                        if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {
                            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                            let mut qseq = record.seq().as_bytes();
                            if record.is_reverse() { revcomp_inplace(&mut qseq) };
                            let index = read_index[&qname];
                            unsafe { // indices retrieved from read_index are supposed to be mutually exclusive
                                read_sequences_ref.write(index, BitSeq::from_utf8(&qseq));
                            }
                        }
                    }
                }
            });
            thread_handles.push(handle);
        }
        thread_handles.into_iter()
            .for_each(|h| { h.join().unwrap(); });
    });

    read_sequences
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
