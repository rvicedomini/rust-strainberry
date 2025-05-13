use std::path::Path;

use anyhow::{bail, Context, Error, Result};
use ahash::AHashMap as HashMap;
use itertools::Itertools;
use rust_htslib::bam::{HeaderView,Read,IndexedReader};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};

use crate::cli::Options;
use crate::bam::BamRecordId;
use crate::seq::SeqDatabase;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

impl MappingType {

    pub fn is_containment(&self) -> bool {
        matches!(self, Self::QueryContained | Self::ReferenceContained)
    }

    pub fn is_dovetail(&self) -> bool {
        matches!(self, Self::DovetailPrefix | Self::DovetailSuffix)
    }
}

pub fn classify_mapping(query_range:(usize,usize,usize), target_range:(usize,usize,usize), overhang:usize, r:f64) -> MappingType {
    use MappingType::*;

    let (b1,e1,l1) = query_range;
    let (b2,e2,l2) = target_range;
    assert!(b1<e1 && e1<=l1 && b2<e2 && e2<=l2);
    let left_overhang = std::cmp::min(b1,b2);
    let right_overhang = std::cmp::min(l1-e1,l2-e2);
    let maplen = std::cmp::max(e1-b1,e2-b2);
    let oh_threshold = std::cmp::min(overhang, ((maplen as f64)*r) as usize);

    if b2 <= b1 && b2 <= oh_threshold && right_overhang > oh_threshold {
        ReferencePrefix
    } else if b1 <= b2 && b1 <= oh_threshold && right_overhang > oh_threshold {
        QueryPrefix
    } else if left_overhang > oh_threshold && l2-e2 <= oh_threshold && l2-e2 <= l1-e1 {
        ReferenceSuffix
    } else if left_overhang > oh_threshold && l1-e1 <= oh_threshold && l1-e1 <= l2-e2 {
        QuerySuffix
    } else if left_overhang > oh_threshold || right_overhang > oh_threshold {
        Internal
    } else if b1 >= b2 && l1-e1 >= l2-e2 {
        ReferenceContained
    } else if b2 >= b1 && l2-e2 >= l1-e1 {
        QueryContained
    } else if b1 <= b2 {
        assert!(l2-e2 <= l1-e1);
        DovetailPrefix
    } else {
        assert!(b2<=b1 && l1-e1 <= l2-e2);
        DovetailSuffix
    }
}


#[derive(Debug)]
pub struct SeqAlignment {
    query_idx: usize,
    query_length: usize,
    query_beg: usize,
    query_end: usize,
    strand: u8,
    target_idx: usize,
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

    pub fn bam_record_id(&self) -> BamRecordId { 
        BamRecordId::new(self.query_idx, self.query_beg, self.query_end)
    }

    pub fn qid(&self) -> usize { self.query_idx }
    pub fn query_index(&self) -> usize { self.query_idx }
    pub fn query_length(&self) -> usize { self.query_length }
    pub fn query_beg(&self) -> usize { self.query_beg }
    pub fn query_end(&self) -> usize { self.query_end }
    
    pub fn strand(&self) -> u8 { self.strand }
    pub fn is_forward(&self) -> bool { self.strand() == b'+' }
    pub fn is_reverse(&self) -> bool { self.strand() == b'-' }

    pub fn tid(&self) -> usize { self.target_idx }
    pub fn target_index(&self) -> usize { self.target_idx }
    pub fn target_length(&self) -> usize { self.target_length }
    pub fn target_beg(&self) -> usize { self.target_beg }
    pub fn target_end(&self) -> usize { self.target_end }

    pub fn mapq(&self) -> u8 { self.mapq }
    pub fn is_secondary(&self) -> bool { self.is_secondary }
    pub fn identity(&self) -> f64 { if self.mapping_length > 0 { (self.matches as f64)/(self.mapping_length as f64) } else { 0.0 } }

    pub fn cigar(&self) -> &[Cigar] { self.cigar.as_slice() }
    pub fn cigar_string(&self) -> String { self.cigar.iter().map(|op| op.to_string()).join("") }

    pub fn from_bam_record(record: &Record, header: &HeaderView, ref_db: &SeqDatabase, read_db: &SeqDatabase) -> SeqAlignment {

        let query_name = std::str::from_utf8(record.qname()).unwrap();
        let query_idx = read_db.get_index(query_name);
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

        let strand = if is_reverse { b'-' } else { b'+' };

        let target_name = std::str::from_utf8(header.tid2name(record.tid() as u32)).unwrap();
        let target_idx = ref_db.get_index(target_name);

        let target_beg = record.reference_start() as usize;
        let target_end = record.reference_end() as usize;
        let target_length = ref_db.sequences[target_idx].len();

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
            _ => panic!("Unexpected type of NM field for bam record {record:?}")
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
            query_idx, query_length, query_beg, query_end,
            strand,
            target_idx, target_length, target_beg, target_end,
            matches: exact_matches, mapping_length, mapq,
            is_secondary,
            cigar
        }
    }

    pub fn split(&self, min_indel:usize) -> Vec<SeqAlignment> {

        let mut alignments = vec![];

        let mut query_beg = if self.is_forward() { self.query_beg() } else { self.query_length() - self.query_end() };
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
                        //query_name: self.query_name().to_string(),
                        query_idx: self.query_index(),
                        query_length: self.query_length(),
                        query_beg: if self.is_forward() { query_beg } else { self.query_length() - query_end },
                        query_end: if self.is_forward() { query_end } else { self.query_length() - query_beg },
                        strand: self.strand(),
                        // target_name: self.target_name().to_string(),
                        target_idx: self.target_index(),
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
                query_idx: self.query_index(),
                query_length: self.query_length(),
                query_beg: if self.is_forward() { query_beg } else { self.query_length() - query_end },
                query_end: if self.is_forward() { query_end } else { self.query_length() - query_beg },
                strand: self.strand(),
                target_idx: self.target_index(),
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
    cigar: &'a [Cigar],
}

impl Iterator for IterAlignedBlocks<'_> {
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


pub fn load_bam_alignments(bam_path: &Path, ref_db: &SeqDatabase, read_db: &SeqDatabase, opts: &Options) -> HashMap<usize,Vec<SeqAlignment>> {

    let mut read_alignments: HashMap<usize, Vec<SeqAlignment>> = HashMap::new();
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

        let seqalign = SeqAlignment::from_bam_record(&record, &bam_header, ref_db, read_db);
        read_alignments.entry(seqalign.query_index())
            .or_default()
            .push(seqalign);
    }

    read_alignments
}


// cigar must not contain clipping/padding operations
pub fn first_match_from(position:usize, mut query_beg:usize, mut target_beg:usize, cigar: &[Cigar]) -> Option<(usize,usize)> {
    let mut cig_i = 0;
    let mut cig_len = 0;
    let mut match_found = false;
    while (target_beg <= position || !match_found) && cig_i < cigar.len() {
        match_found = false;
        let cig_op = cigar[cig_i]; cig_len = cig_op.len() as usize;
        match cig_op {
            Cigar::Match(_)|Cigar::Diff(_)|Cigar::Equal(_) => {
                target_beg += cig_len;
                query_beg += cig_len;
                match_found = true;
            }
            Cigar::Del(_)|Cigar::RefSkip(_) => {
                target_beg += cig_len;
            }
            Cigar::Ins(_) => {
                query_beg += cig_len;
            }
            _ => {
                unreachable!("Unexpected cigar operation")
            }
        }
        cig_i += 1;

    }
    if target_beg > position && match_found {
        let bases_ahead = std::cmp::min(target_beg-position,cig_len);
        assert!(bases_ahead <= query_beg && bases_ahead <= target_beg);
        let query_match_pos = query_beg - bases_ahead;
        let target_match_pos = target_beg - bases_ahead;
        return Some((query_match_pos,target_match_pos))
    }
    None
}

// cigar string must not contain clipped bases at the beginning
pub fn last_match_until(position:usize, mut query_beg:usize, mut target_beg:usize, cigar: &[Cigar]) -> Option<(usize,usize)> {
    let mut cig_i = 0;
    let mut match_found = false;
    let mut query_match_pos = query_beg;
    let mut target_match_pos = target_beg;
    while (target_beg < position || !match_found) && cig_i < cigar.len() {
        match cigar[cig_i] {
            Cigar::Match(cig_len)|Cigar::Diff(cig_len)|Cigar::Equal(cig_len) => {
                let inc = if target_beg + (cig_len as usize) <= position { cig_len as usize } else { position-target_beg };
                target_beg += inc;
                query_beg += inc;
                match_found = true;
                query_match_pos = query_beg;
                target_match_pos = target_beg;
            }
            Cigar::Del(cig_len)|Cigar::RefSkip(cig_len) => {
                target_beg += cig_len as usize;
            }
            Cigar::Ins(cig_len) => {
                query_beg += cig_len as usize;
            }
            _ => {
                unreachable!("Unexpected cigar operation")
            }
        }
        cig_i += 1;

    }
    if match_found {
        return Some((query_match_pos,target_match_pos))
    }
    None
}


#[derive(Debug)]
pub struct PafAlignment {
    pub query_name: String,
    pub query_length: usize,
    pub query_beg: usize,
    pub query_end: usize,
    pub strand: u8,
    pub target_name: String,
    pub target_length: usize,
    pub target_beg: usize,
    pub target_end: usize,
    pub matches: usize, // exact matches
    pub mapping_length: usize,
    pub mapq: u8,
    pub identity: f64,
    pub is_secondary: bool,
    pub cigar: CigarString,
}

impl PafAlignment {

    pub fn map_type(&self, overhang: usize, r: f64) -> MappingType {
        let query_range = if self.strand == b'+' {
            (self.query_beg, self.query_end, self.query_length)
        } else {
            (self.query_length-self.query_end, self.query_length-self.query_beg, self.query_length)
        };
        let target_range = (self.target_beg, self.target_end, self.target_length);
        crate::alignment::classify_mapping(query_range, target_range, overhang, r)
    }
}

impl std::str::FromStr for PafAlignment {

    type Err = Error;
    
    fn from_str(line: &str) -> Result<Self,Self::Err> {

        let cols = line.split('\t').collect_vec();
        if cols.len() < 12 { bail!("cannot parse PAF line (missing fields)") }
        
        let query_name = cols[0].to_string();
        let [query_length, query_beg, query_end] = cols[1..4].iter().map(|s|s.parse::<usize>().unwrap()).collect_array().unwrap();
        
        let strand = match cols[4] {
            "+" => b'+',
            "-" => b'-',
            _ => bail!("unrecognised strand field: {}", cols[4])
        };
        
        let target_name = cols[5].to_string();
        let [target_length, target_beg, target_end] = cols[6..9].iter().map(|s|s.parse::<usize>().unwrap()).collect_array().unwrap();
        
        let matches = cols[9].parse().unwrap();
        let mapping_length = cols[10].parse().unwrap();
        let mapq = cols[11].parse().unwrap();
        let identity = 100.0 * (matches as f64) / (mapping_length as f64);
        
        let tags: HashMap<&str,&str> = cols[12..].iter().filter_map(|s| {
                let [key,_,val] = s.splitn(3,':').collect_array()?;
                Some((key,val))
            }).collect();
        let tp = tags.get("tp").context("missing \"tp\" tag in PAF line")?;
        let is_secondary = *tp == "S";
        let cigar = crate::bam::parse_cigar_bytes(tags.get("cg").unwrap_or(&"").as_bytes());
        
        Ok(PafAlignment {
            query_name,
            query_length,
            query_beg,
            query_end,
            strand,
            target_name,
            target_length,
            target_beg,
            target_end,
            matches,
            mapping_length,
            mapq,
            identity,
            is_secondary,
            cigar
        })
    }
}

impl std::fmt::Display for PafAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let PafAlignment {
            query_name, query_length, query_beg, query_end,
            strand,
            target_name, target_length, target_beg, target_end,
            matches, mapping_length, mapq,
            identity:_, is_secondary:_, cigar:_ } = self;
        write!(f, "{query_name}\t{query_length}\t{query_beg}\t{query_end}\t{}\t{target_name}\t{target_length}\t{target_beg}\t{target_end}\t{matches}\t{mapping_length}\t{mapq}", *strand as char)
    }
}
