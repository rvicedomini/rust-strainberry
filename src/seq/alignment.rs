use std::path::Path;
use std::str::FromStr;

use itertools::Itertools;

use rust_htslib::bam::{HeaderView,Read,IndexedReader};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Record, Cigar, Aux};

use rustc_hash::FxHashMap;

use crate::cli::Options;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse
}

impl FromStr for Strand {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Err(())
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum MappingType {
    DovetailSuffix,
    QueryPrefix,
    QuerySuffix,
    ReferencePrefix,
    ReferenceSuffix,
    Internal,
    QueryContained,
    ReferenceContained,
    DovetailPrefix
}


#[derive(Debug)]
pub struct SeqAlignment {
    query_name: String,
    query_length: usize,
    query_beg: usize,
    query_end: usize,
    strand: Strand,
    target_name: String,
    target_length: usize,
    target_beg: usize,
    target_end: usize,
    matches: usize, // exact matches
    mapping_length: usize,
    mapq: u8,
    is_secondary: bool,
    cigar: Vec<Cigar>,
}


impl SeqAlignment {

    pub fn query_name(&self) -> &str { &self.query_name }
    pub fn query_length(&self) -> usize { self.query_length }
    pub fn query_beg(&self) -> usize { self.query_beg }
    pub fn query_end(&self) -> usize { self.query_end }
    
    pub fn strand(&self) -> Strand { self.strand }

    pub fn target_name(&self) -> &str { &self.target_name }
    pub fn target_length(&self) -> usize { self.target_length }
    pub fn target_beg(&self) -> usize { self.target_beg }
    pub fn target_end(&self) -> usize { self.target_end }

    pub fn mapq(&self) -> u8 { self.mapq }
    pub fn is_secondary(&self) -> bool { self.is_secondary }

    pub fn cigar(&self) -> &Vec<Cigar> { &self.cigar }

    pub fn from_bam_record(record: &Record, header: &HeaderView) -> SeqAlignment {

        let query_name = String::from_utf8_lossy(record.qname()).to_string();
        let query_length = record.seq_len_from_cigar(true);

        let cigar = record.cigar();
        let mut query_beg = match *cigar.first().unwrap() {
            Cigar::HardClip(len) => len as usize,
            Cigar::SoftClip(len) => len as usize,
            _ => 0,
        };
        let mut query_end = query_length - match *cigar.last().unwrap() {
            Cigar::HardClip(len) => len as usize,
            Cigar::SoftClip(len) => len as usize,
            _ => 0,
        };

        let is_reverse = record.is_reverse();
        if is_reverse {
            (query_beg, query_end) = (query_length-query_end, query_length-query_beg);
        }

        let strand = if is_reverse { Strand::Reverse } else { Strand::Forward };

        let tid = record.tid() as u32;
        let target_name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
        let target_beg = record.reference_start() as usize;
        let target_end = record.reference_end() as usize;
        let target_length = header.target_len(tid).unwrap() as usize;

        let (matches,indels) = cigar.iter().fold((0,0), |(matches,indels),op| match *op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => (matches + len as usize, indels),
            Cigar::Del(len) | Cigar::Ins(len) => (matches, indels + len as usize),
            _ => (matches,indels)
        });

        let edit_dist = match record.aux(b"NM").unwrap() {
            Aux::Char(ed) => ed as usize,
            Aux::I8(ed) => ed as usize,
            Aux::U8(ed) => ed as usize,
            Aux::I16(ed) => ed as usize,
            Aux::U16(ed) => ed as usize,
            Aux::I32(ed) => ed as usize,
            Aux::U32(ed) => ed as usize,
            _ => panic!("Cannot read NM field from bam record {record:?}")
        };

        let exact_matches: usize = matches - (edit_dist - indels);
        let mapping_length = matches + indels;

        let mapq = record.mapq();
        let is_secondary = record.is_secondary();

        let cigar = cigar.iter()
            .filter(|op| !matches!(op, Cigar::HardClip(_)|Cigar::SoftClip(_)))
            .cloned()
            .collect_vec();

        SeqAlignment { 
            query_name, query_length, query_beg, query_end,
            strand,
            target_name, target_length, target_beg, target_end,
            matches: exact_matches, mapping_length, mapq,
            is_secondary,
            cigar
        }
    }

    pub fn split(&self, min_indel:usize) -> Vec<SeqAlignment> {

        let mut alignments = vec![];

        let mut query_beg = if self.strand() == Strand::Forward { self.query_beg() } else { self.query_length() - self.query_end() };
        let mut target_beg = self.target_beg();
        let (mut query_end, mut target_end) = (query_beg, self.target_end);
        let (mut matches, mut mapping_length) = (0,0);
        let mut cigar = vec![];

        for cigar_op in self.cigar() {
            match cigar_op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let len = *len as usize;
                    query_end += len;
                    target_end += len;
                    matches += len;
                    mapping_length += len;
                },
                Cigar::Del(len) | Cigar::Ins(len) if (*len as usize) >= min_indel => {
                    alignments.push(SeqAlignment {
                        query_name: self.query_name().to_string(),
                        query_length: self.query_length(),
                        query_beg: if self.strand() == Strand::Forward { query_beg } else { self.query_length() - query_end },
                        query_end: if self.strand() == Strand::Forward { query_end } else { self.query_length() - query_beg },
                        strand: self.strand(),
                        target_name: self.target_name().to_string(),
                        target_length: self.target_length(), target_beg, target_end,
                        matches,
                        mapping_length,
                        mapq: self.mapq(),
                        is_secondary: self.is_secondary(),
                        cigar
                    });
                    query_beg = if let Cigar::Ins(_) = cigar_op { query_end + (*len as usize) } else { query_end };
                    target_beg = if let Cigar::Del(_) = cigar_op { target_end + (*len as usize) } else { target_end };
                    (query_end, target_end) = (query_beg, target_end);
                    (matches, mapping_length) = (0,0);
                    cigar = vec![];
                    continue
                },
                Cigar::Del(len) => {
                    target_end += *len as usize;
                    mapping_length += *len as usize;
                },
                Cigar::Ins(len) => {
                    query_end += *len as usize;
                    mapping_length += *len as usize;
                }
                _ => { continue }
            }
            cigar.push(*cigar_op);
        }

        if !cigar.is_empty() {
            alignments.push(SeqAlignment {
                query_name: self.query_name().to_string(),
                query_length: self.query_length(),
                query_beg: if self.strand() == Strand::Forward { query_beg } else { self.query_length() - query_end },
                query_end: if self.strand() == Strand::Forward { query_end } else { self.query_length() - query_beg },
                strand: self.strand(),
                target_name: self.target_name().to_string(),
                target_length: self.target_length(), target_beg, target_end,
                matches,
                mapping_length,
                mapq: self.mapq(),
                is_secondary: self.is_secondary(),
                cigar
            });
        }

        alignments
    }

    pub fn aligned_blocks(&self) -> IterAlignedBlocks {
        IterAlignedBlocks {
            query_pos: self.query_beg(),
            target_pos: self.target_beg(),
            index: 0,
            cigar: self.cigar(),
        }
    }

}


pub struct IterAlignedBlocks<'a> {
    query_pos: usize,
    target_pos: usize,
    index: usize,
    cigar: &'a Vec<Cigar>,
}

impl<'a> Iterator for IterAlignedBlocks<'a> {
    type Item = [usize;4];
    
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(cigar) = self.cigar.get(self.index) {
            self.index += 1;
            match cigar {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let len = *len as usize;
                    let query_beg = self.query_pos;
                    let query_end = self.query_pos + len;
                    let target_beg = self.target_pos;
                    let target_end = self.target_pos + len;
                    self.query_pos += len;
                    self.target_pos += len;
                    return Some([query_beg, query_end, target_beg, target_end]);
                },
                Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                    self.query_pos += *len as usize;
                },
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    self.target_pos += *len as usize;
                }
                Cigar::Pad(_) => panic!("Cigar should not contain padding operations!")
            }
        }
        None
    }
}


pub fn load_bam_alignments(bam_path: &Path, opts: &Options) -> FxHashMap<String,Vec<SeqAlignment>> {

    let mut read_alignments = FxHashMap::default();
    let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
    let bam_header = bam_reader.header().clone();

    bam_reader.fetch(".").expect("Failed fetching records from BAM");
    for record in bam_reader.rc_records() {
        let record = record.expect("Failed processing BAM file");
        if record.mapq() < opts.min_mapq
            || record.is_unmapped()
            || record.is_secondary()
            || record.is_quality_check_failed()
            || record.is_duplicate()
        {
            continue
        }

        let seqalign = SeqAlignment::from_bam_record(&record, &bam_header);
        for seqalign in seqalign.split(opts.min_indel) {
            if !read_alignments.contains_key(seqalign.query_name()) {
                read_alignments.insert(seqalign.query_name().to_string(), vec![]);
            }
            read_alignments.get_mut(seqalign.query_name()).unwrap().push(seqalign);
        }
    }

    read_alignments

}