pub mod alignment;
pub mod read;

use std::fmt;
use std::ops::Range;
// use std::path::Path;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::utils::BamRecordId;


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
        write!(f, "{}: {} ({}:{}..={})", self.record_id().0, seq_string, self.target_id(), self.positions().first().unwrap(), self.positions().last().unwrap())
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
    pub fn target_id(&self) -> usize { self.target_id }
    pub fn positions(&self) -> &Vec<usize> { &self.positions }
    pub fn nucleotides(&self) -> &Vec<u8> { &self.nucleotides }

    pub fn len(&self) -> usize { self.positions.len() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn range(&self) -> Range<usize> {
        let beg = *self.positions.first().unwrap();
        let end = *self.positions.last().unwrap();
        beg..end
    }

    pub fn seq_interval(&self) -> SeqInterval {
        SeqInterval { 
            tid: self.target_id,
            beg: *self.positions.first().unwrap(),
            end: *self.positions.last().unwrap() + 1,
        }
    }

    pub fn push(&mut self, pos:usize, nuc:u8) {
        self.positions.push(pos);
        self.nucleotides.push(nuc);
    }

    pub fn from_bam_record(record: &bam::record::Record, variant_positions: &[usize], mut var_idx: usize) -> Option<SuccinctSeq> {
        
        if var_idx >= variant_positions.len() {
            return None
        }

        let record_id = crate::utils::bam_record_id(record);
        let target_id = record.tid() as usize;
        let mut sseq = SuccinctSeq::build(record_id, target_id);

        let mut target_pos = record.pos() as usize;
        let mut query_pos = 0;
        let query_seq = record.seq();

        let mut var_pos = variant_positions[var_idx];
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
                        if var_idx >= variant_positions.len() {
                            break;
                        }
                        var_pos = variant_positions[var_idx];
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
                        if var_idx >= variant_positions.len(){
                            break;
                        }
                        var_pos = variant_positions[var_idx];
                        target_pos += 1;
                        oplen -= 1;
                    }
                },
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }

            if var_idx >= variant_positions.len() {
                break;
            }
        }

        Some(sseq)
    }
}



// pub fn build_succinct_sequences(bam_path: &Path, target_names: &Vec<String>, target_intervals: &'a Vec<SeqInterval>, read_sequences: &'a FxHashMap<String, Vec<u8>>, output_dir: &Path, opts: &'a Options) {
//     todo!()
// }


// def build_succinct_reads(bamfile, reference_id, reference_positions, min_mapq):
//     succinct_reads = dict()
//     with pysam.AlignmentFile(bamfile,'rb') as bam_handle:
//         for segment in bam_handle.fetch(contig=reference_id):
//             if segment.is_unmapped or segment.is_secondary or segment.is_duplicate or segment.is_qcfail or segment.mapping_quality < min_mapq: 
//                 continue
//             pos_idx = bisect.bisect_left(reference_positions,segment.reference_start)
//             sread = get_sread_from_segment(segment,reference_positions,pos_idx)
//             if sread and len(sread.positions) > 0:
//                 succinct_reads[sread.id] = sread
//     return succinct_reads


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