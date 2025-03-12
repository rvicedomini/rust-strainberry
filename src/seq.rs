pub mod alignment;
pub mod read;

use std::fmt;
use std::path::Path;
use std::sync::mpsc;
use std::thread;

use rust_htslib::bam::{Read,IndexedReader};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, Record};

use crate::cli::Options;
use crate::utils::{self,BamRecordId};
use crate::variant::{Var,VarDict};


#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
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
    target_id: usize,
    positions: Vec<usize>,
    nucleotides: Vec<u8>
}


impl fmt::Display for SuccinctSeq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq_string = String::from_utf8_lossy(self.nucleotides()).to_string();
        write!(f, "{}: {} ({}:{}..={})", self.record_id().0, seq_string, self.tid(), self.positions().first().unwrap(), self.positions().last().unwrap())
    }
}


impl SuccinctSeq {

    pub fn build(record_id:BamRecordId, target_id:usize) -> Self {
        Self { 
            record_id,
            target_id,
            positions: vec![],
            nucleotides: vec![] 
        }
    }

    pub fn record_id(&self) -> BamRecordId { self.record_id.clone() }
    pub fn tid(&self) -> usize { self.target_id }
    pub fn positions(&self) -> &Vec<usize> { &self.positions }
    pub fn nucleotides(&self) -> &Vec<u8> { &self.nucleotides }

    pub fn len(&self) -> usize { self.positions.len() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn beg(&self) -> usize { *self.positions.first().unwrap() }
    pub fn end(&self) -> usize { *self.positions.last().unwrap() + 1 }
    pub fn region(&self) -> SeqInterval {
        SeqInterval { tid: self.tid(), beg: self.beg(), end: self.end() }
    }

    pub fn push(&mut self, pos:usize, nuc:u8) {
        self.positions.push(pos);
        self.nucleotides.push(nuc);
    }

    pub fn from_bam_record(record: &Record, target_variants: &[Var], mut var_idx: usize) -> Option<SuccinctSeq> {
        
        if var_idx >= target_variants.len() {
            return None
        }

        let record_id = crate::utils::bam_record_id(record);
        let target_id = record.tid() as usize;
        let mut sseq = SuccinctSeq::build(record_id, target_id);

        let mut target_pos = record.pos() as usize;
        let mut query_pos = 0;
        let query_seq = record.seq();

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
                        sseq.push(target_pos, query_seq[query_pos]);
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
                        sseq.push(target_pos, b'-');
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



pub fn build_succinct_sequences(bam_path: &Path, variants: &VarDict, opts: &Options) -> Vec<SuccinctSeq> {
    
    let mut target_intervals = utils::bam_target_intervals(bam_path);
    target_intervals.sort_unstable_by_key(|siv| siv.end - siv.beg);
    
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
                            if let Some(sseq) = SuccinctSeq::from_bam_record(&record, target_variants, var_idx) {
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


// # TODO: improve it, considering a sorted list of positions in input
// def get_sread_from_pafalignment(a:ReadAlignment, query_sequence:str, variant_positions:set):
//     shtig = SuccinctRead(a.uid(),a.reference)
//     ref_pos = a.reference_start
//     query_pos = a.query_start if a.strand == '+' else a.query_length-a.query_end
//     for cigop, oplen in a.cigar:
//         assert(oplen > 0)
//         if cigop in [CigarOp.MATCH,CigarOp.MATCH_EQ,CigarOp.MISMATCH]:
//             while oplen > 0:
//                 if ref_pos in variant_positions:
//                     shtig.positions.append(ref_pos)
//                     nuc = query_sequence[query_pos] if a.strand == '+' else query_sequence[a.query_length-query_pos-1].translate(RC_TABLE)
//                     shtig.nucleotides.append(nuc)
//                 ref_pos += 1
//                 query_pos += 1
//                 oplen -= 1
//         elif cigop in [CigarOp.DELETION,CigarOp.REF_SKIP]:
//             while oplen > 0:
//                 if ref_pos in variant_positions:
//                     shtig.positions.append(ref_pos)
//                     shtig.nucleotides.append('-')
//                 ref_pos += 1
//                 oplen -= 1
//         elif cigop == CigarOp.INSERTION:
//             query_pos += oplen
//         else:
//             assert(False)
//     return shtig