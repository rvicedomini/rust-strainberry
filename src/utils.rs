use std::fs::File;
use std::io::{BufRead,BufReader};
use std::mem::MaybeUninit;
use std::path::Path;

use flate2::read::MultiGzDecoder;
use needletail::Sequence;
use rustc_hash::FxHashMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::CigarString;

pub mod seq;

pub fn get_maxrss() -> f64 {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as f64 / (1024.0 * 1024.0)
}


pub fn get_file_reader(path: &Path) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap();
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    }
}


pub fn chrom2tid(bam_path:&Path) -> FxHashMap<String,usize> {
    let bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    bam_reader.header()
        .target_names().iter()
        .enumerate()
        .map(|(tid,&name)| (String::from_utf8_lossy(name).to_string(),tid as usize))
        .collect()
}

pub fn load_sequences(fasta_path: &Path, bam_path: &Path) -> Vec<Vec<u8>> {
    let bam_reader = bam::Reader::from_path(bam_path).unwrap();
    let header_view = bam_reader.header();

    let mut target_sequences = vec![Vec::new(); header_view.target_count() as usize];

    let mut fasta_reader = needletail::parse_fastx_file(fasta_path).expect("Cannot open fasta file");
    while let Some(record) = fasta_reader.next() {
        let record = record.unwrap();
        let tid = header_view.tid(record.id()).unwrap() as usize;
        target_sequences[tid] = record.normalize(false).to_vec();
    }

    target_sequences
}


pub fn parse_cigar_bytes(cigar: &[u8]) -> CigarString {
    CigarString::try_from(cigar)
        .expect("Unable to parse cigar string.")
}


pub fn estimate_lookback(bam_path: &Path, n: usize) -> Option<usize> {
    let mut bam_reader = bam::Reader::from_path(bam_path).unwrap();

    let mut read_lengths = Vec::with_capacity(n);
    
    let mut record = bam::Record::new();
    while read_lengths.len() < n && bam_reader.read(&mut record).is_some() {
        
        if record.is_unmapped() 
            || record.is_secondary() 
            || record.is_quality_check_failed() 
            || record.is_duplicate() 
            || record.is_supplementary()
        {
            continue;
        }

        read_lengths.push(record.seq_len());
    }

    if read_lengths.len() == 0 {
        return None
    }

    read_lengths.sort_unstable();
    
    let lookup = read_lengths[read_lengths.len()/10];
    Some(lookup)
}