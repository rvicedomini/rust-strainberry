use std::fmt;
use std::path::Path;
use std::sync::mpsc;
use std::thread;

use ahash::AHashMap as HashMap;
use anyhow::{bail,Result};
// use bitseq::BitSeq;
use itertools::Itertools;
use needletail::Sequence;
use rust_htslib::bam::{Read,IndexedReader};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, Record};

use crate::bitseq::BitSeq;
use crate::cli::Options;
use crate::bam::BamRecordId;
use crate::variant::{Var,VarDict};


// flip strand as u8 character
// b'+': 00101|01|1
// b'-': 00101|10|1
//    6: 00000|11|0
#[inline(always)]
pub fn flip_strand(strand:u8) -> u8 {
    assert!(strand == b'+' || strand == b'-');
    strand ^ 6
}

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


#[derive(Default)]
pub struct SeqDatabase {
    pub index: HashMap<String,usize>,
    pub sequences: Vec<BitSeq>,
    pub names: Vec<String>,
    pub nb_bases: usize,
}

impl SeqDatabase {

    pub fn build(path: &Path, index_names: bool) -> Result<SeqDatabase> {

        let mut reader = needletail::parse_fastx_file(path)?;
    
        let mut index = HashMap::new();
        let mut sequences = Vec::new();
        let mut names = Vec::new();
    
        let mut nb_bases = 0;
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let name = std::str::from_utf8(record.id().split(|b| b.is_ascii_whitespace()).next().unwrap()).unwrap();
            if index.insert(name.to_string(), index.len()).is_some() {
                bail!("duplicate reference identifier: {name}");
            }
            let sequence = BitSeq::from_utf8(record.normalize(false).as_ref());
            nb_bases += sequence.len();
            sequences.push(sequence);
            if index_names {
                names.push(name.to_string());
            }
        }
        assert!(index.len() == sequences.len());
        Ok(SeqDatabase { index, sequences, names, nb_bases })
    }

    pub fn size(&self) -> usize {
        self.sequences.len()
    }

    pub fn get_index(&self, name: &str) -> usize {
        self.index[name]
    }

    pub fn compute_nx(&self, x:usize) -> Option<usize> {
        assert!((10..=90).contains(&x));

        let mut lengths = Vec::with_capacity(self.sequences.len());
        self.sequences.iter().for_each(|seq| lengths.push(seq.len()));
        lengths.sort_unstable_by(|a,b| b.cmp(a));

        let mut cumulative_len = 0;
        for len in lengths {
            cumulative_len += len;
            if 100 * cumulative_len >= x * self.nb_bases {
                return Some(len)
            }
        }
        None
    }
}


#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeqInterval {
    pub tid: usize,
    pub beg: usize,
    pub end: usize,
}

impl fmt::Display for SeqInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.tid, self.beg, self.end)
    }
}

impl SeqInterval {
    pub fn length(&self) -> usize {
        self.end - self.beg
    }

    pub fn empty(&self) -> bool {
        self.end == self.beg
    }
}


pub struct SuccinctSeq {
    record_id: BamRecordId,
    ref_idx: usize,
    positions: Vec<usize>,
    nucleotides: Vec<u8>,
    qualities: Vec<u8>,
}


impl fmt::Display for SuccinctSeq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq_string = String::from_utf8_lossy(self.nucleotides()).to_string();
        write!(f, "{}: {} ({}:{}..={})", self.record_id().index, seq_string, self.tid(), self.positions().first().unwrap(), self.positions().last().unwrap())
    }
}


impl SuccinctSeq {

    pub fn new(record_id:BamRecordId, ref_idx:usize) -> Self {
        Self { 
            record_id,
            ref_idx,
            positions: Vec::new(),
            nucleotides: Vec::new(),
            qualities: Vec::new(),
        }
    }

    pub fn record_id(&self) -> BamRecordId { self.record_id }
    pub fn tid(&self) -> usize { self.ref_idx }
    pub fn positions(&self) -> &[usize] { &self.positions }
    pub fn nucleotides(&self) -> &[u8] { &self.nucleotides }
    pub fn qualities(&self) -> &[u8] { &self.qualities }

    pub fn len(&self) -> usize { self.positions.len() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn beg(&self) -> usize { *self.positions.first().unwrap() }
    pub fn end(&self) -> usize { *self.positions.last().unwrap() + 1 }
    pub fn region(&self) -> SeqInterval {
        SeqInterval { tid: self.tid(), beg: self.beg(), end: self.end() }
    }

    pub fn push(&mut self, pos:usize, nuc:u8, qual:u8) {
        self.positions.push(pos);
        self.nucleotides.push(nuc);
        self.qualities.push(qual);
    }

    pub fn from_bam_record(record: &Record, target_variants: &[Var], mut var_idx: usize, read_db: &SeqDatabase) -> Option<SuccinctSeq> {
        
        if var_idx >= target_variants.len() {
            return None
        }

        let record_id = BamRecordId::from_record(record, read_db);
        let target_id = record.tid() as usize;
        let mut sseq = SuccinctSeq::new(record_id, target_id);

        let mut target_pos = record.pos() as usize;
        let mut query_pos = 0;
        let query_seq = record.seq();
        let query_qual = record.qual();

        let mut var_pos = target_variants[var_idx].pos;
        for cig in record.cigar().iter() {
            match cig {
                Cigar::SoftClip(len) | Cigar::Ins(len) => {
                    query_pos += *len as usize;
                    continue;
                },
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let mut oplen = *len as usize;
                    while oplen > 0 {
                        assert!(var_pos >= target_pos);
                        if var_pos >= target_pos + oplen {
                            target_pos += oplen;
                            query_pos += oplen;
                            break
                        }
                        // jump to variant position
                        let dist = var_pos - target_pos;
                        target_pos += dist;
                        query_pos += dist;
                        oplen -= dist;
                        sseq.push(target_pos, query_seq[query_pos], query_qual[query_pos]);
                        // consume match on both reference/query and retrieve next variant position
                        var_idx += 1;
                        if var_idx >= target_variants.len() {
                            break;
                        }
                        var_pos = target_variants[var_idx].pos;
                        target_pos += 1;
                        query_pos += 1;
                        oplen -= 1;
                    }
                },
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    let mut oplen = *len as usize;
                    while oplen > 0 {
                        assert!(var_pos >= target_pos);
                        if var_pos >= target_pos + oplen {
                            target_pos += oplen;
                            break;
                        }
                        // jump to variant position
                        let dist = var_pos - target_pos;
                        target_pos += dist;
                        oplen -= dist;
                        sseq.push(target_pos, b'-', 0);
                        // consume deletion on reference and retrieve next variant position
                        var_idx += 1;
                        if var_idx >= target_variants.len(){
                            break;
                        }
                        var_pos = target_variants[var_idx].pos;
                        target_pos += 1;
                        oplen -= 1;
                    }
                },
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }

            if var_idx >= target_variants.len() {
                break;
            }
        }

        if !sseq.is_empty() {
            Some(sseq)
        } else {
            None
        }
    }
}


pub fn build_succinct_sequences(bam_path: &Path, ref_db: &SeqDatabase, read_db: &SeqDatabase, variants: &VarDict, opts: &Options) -> Vec<SuccinctSeq> {

    let index_to_tid = crate::bam::build_target_index(bam_path, ref_db);
    let target_intervals = ref_db.sequences.iter().enumerate()
        .map(|(ref_idx,seq)| SeqInterval{ tid:index_to_tid[ref_idx], beg:0, end:seq.len() })
        .sorted_unstable_by_key(|iv| -((iv.end-iv.beg) as isize))
        .collect_vec();
    
    let (tx, rx) = mpsc::channel();
    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let target_intervals_ref = &target_intervals;
            scope.spawn(move || {
                let mut succinct_sequences = vec![];
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for &siv in target_intervals_ref.iter().skip(thread_id).step_by(opts.nb_threads) {
                    if let Some(target_variants) = variants.get(&siv.tid) {

                        bam_reader.fetch(siv.tid as u32).unwrap();
                        for record in bam_reader.rc_records().map(|x| x.expect("Failure parsing BAM file")) {
                            if record.mapq() < opts.min_mapq
                                || record.is_unmapped()
                                || record.is_secondary()
                                || record.is_quality_check_failed()
                                || record.is_duplicate()
                            {
                                continue
                            }
                            let var_idx = target_variants.partition_point(|var| var.pos < record.reference_start() as usize);
                            if let Some(sseq) = SuccinctSeq::from_bam_record(&record, target_variants, var_idx, read_db) {
                                succinct_sequences.push(sseq);
                            }
                        }
                    }
                }
                sender.send(succinct_sequences).unwrap();
            });
        }
    });

    (0..opts.nb_threads)
        .flat_map(|_| rx.recv().unwrap().into_iter())
        .collect()
}
