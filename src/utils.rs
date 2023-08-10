use std::fs::File;
use std::io::{BufRead,BufReader,Write,BufWriter};
use std::mem::MaybeUninit;
use std::path::Path;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use itertools::Itertools;
use needletail::Sequence;
use rustc_hash::FxHashMap;

use rust_htslib::bam::{self, HeaderView};
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar,CigarString};

use crate::seq::SeqInterval;


pub fn get_maxrss() -> f64 {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as f64 / (1024.0 * 1024.0)
}


pub fn get_file_reader(path: &Path) -> Box<dyn BufRead> {
    let file = match File::open(path) {
        Err(err) => panic!("Could not open file \"{}\": {}", path.display(), err),
        Ok(file) => file,
    };
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    }
}


pub fn get_file_writer(path: &Path) -> Box<dyn Write> {
    let file = match File::create(&path) {
        Err(err) => panic!("Could not open/create file \"{}\": {}", path.display(), err),
        Ok(file) => file,
    };
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(
            BufWriter::new(GzEncoder::new(file, Compression::default()))
        ),
        _ => Box::new(BufWriter::new(file)),
    }
}


pub fn bam_target_names(bam_path:&Path) -> Vec<String> {
    let bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    bam_reader
        .header()
        .target_names()
        .iter()
        .map(|&name| String::from_utf8_lossy(name).to_string())
        .map(|x| x)
        .collect_vec()
}


pub fn bam_target_intervals(bam_path:&Path) -> Vec<SeqInterval> {
    
    let bam_reader = bam::Reader::from_path(bam_path).unwrap();
    let bam_header = bam_reader.header();

    (0..bam_header.target_count())
        .map(|tid| SeqInterval {
            tid: tid as usize,
            beg: 0,
            end: bam_header.target_len(tid).unwrap() as usize
        }).collect_vec()
}


pub fn chrom2tid(bam_header:&HeaderView) -> FxHashMap<String,usize> {
    bam_header
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
        // eprintln!("{tid} -> {}", String::from_utf8_lossy(record.id()));
        target_sequences[tid] = record.normalize(false).to_vec();
    }

    target_sequences
}


pub fn parse_cigar_bytes(cigar: &[u8]) -> CigarString {
    CigarString::try_from(cigar)
        .expect("Unable to parse cigar string.")
}


pub fn intervals_from_cigar(cigarstring: &CigarString, pos:usize) -> [usize;4] {

    let (target_beg, mut target_end) = (pos,pos);
    let (mut query_beg, mut query_end) = (0,0);
    for (i,&cigar) in cigarstring.iter().enumerate() {
        match cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                target_end += len as usize;
                query_end += len as usize;
            },
            Cigar::Ins(len) => {
                query_end += len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                target_end += len as usize;
            }
            Cigar::SoftClip(len) | Cigar::HardClip(len) if i == 0 => {
                query_beg += len as usize;
                query_end += len as usize;
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                break;
            }
            Cigar::Pad(_) => panic!("Cigar should not contain padding operations!")
        }
    }

    [ query_beg, query_end, target_beg, target_end ]
}


pub fn seq_length_from_cigar(cigarstring: &CigarString, include_hard_clip: bool) -> usize {
    let mut length = 0;
    for &cigar in cigarstring.iter() {
        match cigar {
            Cigar::Match(len) | Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                length += len as usize;
            }
            Cigar::HardClip(len) if include_hard_clip => {
                length += len as usize;
            }
            _ => {}
        }
    }
    length
}


#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct BamRecordId(pub String, pub usize, pub usize);

pub fn bam_record_id(record: &bam::record::Record) -> BamRecordId {
    let cigar = record.cigar();
    
    let query_name = String::from_utf8(record.qname().to_vec()).unwrap();
    let query_length = record.seq_len_from_cigar(true);

    let mut query_start = match *cigar.first().unwrap() {
        Cigar::HardClip(len) => len as usize,
        Cigar::SoftClip(len) => len as usize,
        _ => 0,
    };

    let mut query_end = query_length - match *cigar.last().unwrap() {
        Cigar::HardClip(len) => len as usize,
        Cigar::SoftClip(len) => len as usize,
        _ => 0,
    };

    if record.is_reverse() {
        (query_start, query_end) = (query_length-query_end, query_length-query_start);
    }

    BamRecordId(query_name, query_start, query_end)
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


// def hamming_distance(a_seq:str, b_seq:str) -> int:
//     assert(len(a_seq)==len(b_seq))
//     ndiff = sum(a!=b for a,b in zip(a_seq,b_seq))
//     return ndiff