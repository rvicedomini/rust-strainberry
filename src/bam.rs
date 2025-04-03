use std::path::Path;

use anyhow::{anyhow, Result};
use ahash::AHashMap as HashMap;
use itertools::Itertools;
use needletail::Sequence;
use rust_htslib::bam::{self, HeaderView};
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar,CigarString};
use rust_htslib::htslib::{BAM_FSECONDARY,BAM_FSUPPLEMENTARY};
use tinyvec::{tiny_vec,TinyVec};

use crate::seq::bitseq::BitSeq;
use crate::seq::SeqInterval;


pub fn bam_to_fasta(bam_path:&Path, out_path:&Path) -> Result<()> {

    let mut bam = bam::IndexedReader::from_path(bam_path)?;
    bam.fetch(".").unwrap();

    let primary_iter = bam.rc_records()
        .map(|x| x.expect("Error parsing BAM file"))
        .filter(|rec| rec.flags() & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY) as u16 == 0);

    let mut writer = crate::utils::get_file_writer(out_path);
    for rec in primary_iter {
        writer.write_all(format!(">{}\n", std::str::from_utf8(rec.qname())?).as_bytes())?;
        writer.write_all(&rec.seq().as_bytes())?;
        writer.write_all(b"\n")?;
    }

    Ok(())
}


pub fn bam_target_names(bam_path:&Path) -> (Vec<String>,HashMap<String,usize>) {
    let bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    let target_names = bam_reader.header()
        .target_names()
        .iter()
        .map(|&name| String::from_utf8_lossy(name).to_string())
        .collect_vec();
    let name_to_index = target_names.iter()
        .enumerate()
        .map(|(idx,name)| (name.clone(),idx))
        .collect::<HashMap<String,usize>>();
    (target_names,name_to_index)
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


pub fn chrom2tid(bam_header:&HeaderView) -> HashMap<String,usize> {
    bam_header
        .target_names().iter()
        .enumerate()
        .map(|(target_index, name)| (
            std::str::from_utf8(name).unwrap().to_string(),
            target_index
        )).collect()
}


pub fn load_sequences(fasta_path: &Path, bam_path: &Path) -> Vec<BitSeq> {
    let bam_reader = bam::Reader::from_path(bam_path).unwrap();
    let header_view = bam_reader.header();

    let mut target_sequences = Vec::new();
    target_sequences.resize_with(header_view.target_count() as usize, BitSeq::default);

    let mut fasta_reader = needletail::parse_fastx_file(fasta_path).expect("Cannot open fasta file");
    while let Some(record) = fasta_reader.next() {
        let record = record.unwrap();
        let tid = header_view.tid(record.id()).unwrap() as usize;
        // eprintln!("{tid} -> {}", String::from_utf8_lossy(record.id()));
        target_sequences[tid] = BitSeq::from_utf8(record.normalize(false).as_ref());
    }

    target_sequences
}


pub fn parse_cigar_bytes(cigar: &[u8]) -> CigarString {
    CigarString::try_from(cigar)
        .expect("Unable to parse cigar string.")
}


pub fn intervals_from_cigar(cigarstring: &CigarString, target_beg:usize, mut query_beg:usize) -> [usize;4] {
    let mut target_end = target_beg;
    let mut query_end = query_beg;
    for (i, &cigar) in cigarstring.iter().enumerate() {
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

    if read_lengths.is_empty() {
        return None
    }

    read_lengths.sort_unstable();
    
    let lookup = read_lengths[read_lengths.len()/10];
    Some(lookup)
}


#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct BamRecordId {
    pub index:usize,
    pub beg:usize,
    pub end:usize
}

impl BamRecordId {

    pub fn new(index:usize, beg:usize, end:usize) -> Self {
        Self { index, beg, end }
    }

    pub fn from_record(record: &bam::record::Record, read_index: &HashMap<String,usize>) -> Self {
        let cigar = record.cigar();
    
        let query_name = std::str::from_utf8(record.qname()).unwrap();
        let query_index = read_index[query_name];
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

        // start/end positions should be relative to the forward strand of the read
        if record.is_reverse() {
            (query_start, query_end) = (query_length-query_end, query_length-query_start);
        }

        Self { index:query_index, beg:query_start, end:query_end }
    }
    
}


// mimics IndexedReads of pysam python library
// https://github.com/pysam-developers/pysam/blob/a7871992153506629c435d62e6150d62bb6537d2/pysam/libcalignmentfile.pyx#L2876
pub struct BamReadIndex {
    bam: bam::IndexedReader,
    qname_to_voffsets: HashMap<String,TinyVec<[i64;4]>>,
    fetched: TinyVec<[i64;4]>,
    idx: usize,
}

impl BamReadIndex {

    pub fn build_from(bam_path:&Path) -> BamReadIndex {
        
        let mut bam  = bam::IndexedReader::from_path(bam_path).expect("Error opening BAM file.");
        let mut qname_to_voffsets: HashMap<String,TinyVec<[i64;4]>> = HashMap::new();

        let mut voffset = bam.tell();
        let mut record = bam::Record::new();
        while let Some(result) = bam.read(&mut record) {
            result.expect("BAM parsing failed");
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            qname_to_voffsets.entry(qname).or_default().push(voffset);
            voffset = bam.tell();
        }

        BamReadIndex { bam, qname_to_voffsets, fetched:tiny_vec![], idx:0 }
    }

    pub fn fetch(&mut self, qname:&str) {
        self.idx = 0;
        self.fetched.clear();
        if let Some(voffsets) = self.qname_to_voffsets.get(qname) {
            self.fetched.extend_from_slice(voffsets);
        }
    }

    pub fn read(&mut self, record:&mut bam::Record) -> Result<()> {
        if let Some(offset) = self.fetched.get(self.idx) {
            self.bam.seek(*offset).unwrap();
            if let Some(Ok(())) = self.bam.read(record) {
                self.idx += 1;
                return Ok(());
            };
        }
        self.idx = self.fetched.len();
        Err(anyhow!("no more records"))
    }

    pub fn read_primary(&mut self, qname:&str, record:&mut bam::Record) -> Result<()> {
        self.fetch(qname);
        while self.read(record).is_ok() {
            if record.flags() & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY) as u16 == 0 {
                return Ok(());
            }
        }
        Err(anyhow!("no primary record found"))
    }
}
