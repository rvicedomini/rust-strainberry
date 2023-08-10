pub mod alignment;
pub mod read;

use std::fmt;
use std::ops::Range;

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

pub struct SuccinctSeq {
    name: String,
    tid: usize,
    positions: Vec<usize>,
    nucleotides: Vec<u8>
}

impl SuccinctSeq {

    pub fn build(name:&str, tid:usize) -> Self {
        Self { 
            name: name.to_string(), 
            tid,
            positions: vec![],
            nucleotides: vec![] 
        }
    }

    pub fn name(&self) -> &str { self.name.as_str() }
    pub fn tid(&self) -> usize { self.tid }
    pub fn positions(&self) -> &Vec<usize> { &self.positions }
    pub fn nucleotides(&self) -> &Vec<u8> { &self.nucleotides }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn range(&self) -> Range<usize> {
        let beg = *self.positions.first().unwrap();
        let end = *self.positions.last().unwrap();
        beg..end
    }

    pub fn seq_interval(&self) -> SeqInterval {
        SeqInterval { 
            tid: self.tid,
            beg: *self.positions.first().unwrap(),
            end: *self.positions.last().unwrap() + 1,
        }
    }

    pub fn push(&mut self, pos:usize, nuc:u8) {
        self.positions.push(pos);
        self.nucleotides.push(nuc);
    }
}

impl fmt::Display for SuccinctSeq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq_string = String::from_utf8_lossy(self.nucleotides()).to_string();
        write!(f, "{}: {} ({}:{}..={})", self.name(), seq_string, self.tid(), self.positions().first().unwrap(), self.positions().last().unwrap())
    }
}


// def get_sread_from_segment(segment:pysam.AlignedSegment, variant_positions:list, pos_idx=0) -> SuccinctRead:
//     if pos_idx >= len(variant_positions):
//         return None
//     sread = SuccinctRead(segment_uid(segment), segment.reference_name)
//     ref_pos = segment.reference_start
//     seg_pos = segment.query_alignment_start
//     var_idx = pos_idx
//     var_pos = variant_positions[var_idx]
//     for cigop, oplen in segment.cigartuples:
//         if var_idx >= len(variant_positions): 
//             break
//         if cigop in [CigarOp.MATCH,CigarOp.MATCH_EQ,CigarOp.MISMATCH]:
//             while oplen > 0:
//                 assert(var_pos >= ref_pos)
//                 if var_pos-ref_pos >= oplen:
//                     ref_pos += oplen
//                     seg_pos += oplen
//                     oplen = 0
//                     continue
//                 # jump to variant position
//                 dist = var_pos-ref_pos
//                 ref_pos += dist
//                 seg_pos += dist
//                 oplen -= dist
//                 sread.positions.append(ref_pos)
//                 sread.nucleotides.append(segment.query_sequence[seg_pos])
//                 # consume match on both reference/query and retrieve next variant position
//                 var_idx += 1
//                 if var_idx >= len(variant_positions): break
//                 var_pos = variant_positions[var_idx]
//                 ref_pos += 1
//                 seg_pos += 1
//                 oplen -= 1
//         elif cigop in [CigarOp.DELETION,CigarOp.REF_SKIP]:
//             while oplen > 0:
//                 assert(var_pos >= ref_pos)
//                 if var_pos-ref_pos >= oplen:
//                     ref_pos += oplen
//                     oplen = 0
//                     continue
//                 # jump to variant position
//                 dist = var_pos-ref_pos
//                 ref_pos += dist
//                 oplen -= dist
//                 sread.positions.append(ref_pos)
//                 sread.nucleotides.append('-')
//                 # consume deletion on reference and retrieve next variant position
//                 var_idx += 1
//                 if var_idx >= len(variant_positions): break
//                 var_pos = variant_positions[var_idx]
//                 ref_pos += 1
//                 oplen -= 1
//         elif cigop == CigarOp.INSERTION:
//             seg_pos += oplen
//     return sread


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