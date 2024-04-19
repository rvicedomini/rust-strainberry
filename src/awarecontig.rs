use std::fmt;

use itertools::Itertools;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::phase::haplotype::{HaplotypeId,Haplotype};
use crate::seq::SeqInterval;
use crate::seq::alignment::{MappingType, SeqAlignment, Strand};
use crate::utils::BamRecordId;


#[derive(Debug, Clone)]
pub struct AwareContig {
    interval: SeqInterval,
    is_separated: bool,
    hid: Option<usize>,
    sequence: Vec<u8>,
    aligned_bases: usize,
}

impl fmt::Display for AwareContig {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "AwareContig(iv={}, is_separated={}, ht_id={:?})", self.interval, self.is_separated, self.hid)
    }
}

impl AwareContig {

    fn tid(&self) -> usize { self.interval.tid }
    fn beg(&self) -> usize { self.interval.beg }
    fn end(&self) -> usize { self.interval.end }
    fn length(&self) -> usize { self.interval.end - self.interval.beg }

    fn is_phased(&self) -> bool { self.hid.is_some() }
    fn hid(&self) -> Option<usize> { self.hid }
    fn haplotype_id(&self) -> Option<HaplotypeId> {
        Some(HaplotypeId{
            tid: self.interval.tid,
            beg: self.interval.beg,
            end: self.interval.end,
            hid: self.hid?
        })
    }

    fn interval(&self) -> SeqInterval { self.interval }

    fn phaseset_id(&self) -> SeqInterval { self.interval() }

    fn region(&self) -> (usize,usize) { (self.beg(), self.end()) }

    fn depth(&self) -> f64 {
        let contig_length = self.length();
        if contig_length > 0 {
            (self.aligned_bases as f64) / (contig_length as f64)
        } else {
            0.0
        }
    }

}

#[derive(Debug, Clone)]
pub struct AwareAlignment {
    query_name: String,
    query_length: usize,
    query_beg: usize,
    query_end: usize,
    strand: Strand,
    aware_id: usize,
    target_name: String,
    target_beg: usize,
    target_end: usize,
    mapq: u8,
    identity: f64,
    query_aln_beg: usize,
    query_aln_end: usize,
    mapping_type: MappingType
}


fn retrieve_unphased_intervals(target_intervals:&Vec<SeqInterval>, haplotypes:&FxHashMap<HaplotypeId,Haplotype>, min_length:usize) -> Vec<SeqInterval> {

    let mut unphased_intervals: Vec<SeqInterval> = vec![];
    let phased_intervals = haplotypes.values().map(|ht| ht.region()).unique().sorted_unstable().collect_vec();

    for iv in target_intervals {

        let first = phased_intervals.partition_point(|x| x < iv);
        
        let mut beg = iv.beg;
        for &piv in phased_intervals[first..].iter().take_while(|piv| piv.tid == iv.tid && piv.end <= iv.end) {
            let uiv = SeqInterval{ tid:iv.tid, beg, end:piv.beg };
            if !uiv.empty() && uiv.length() >= min_length {
                unphased_intervals.push(uiv);
            }
            beg = piv.end;
        }

        let uiv = SeqInterval{ tid:iv.tid, beg, end:iv.end };
        if !uiv.empty() && uiv.length() >= min_length {
            unphased_intervals.push(uiv);
        }
    }

    unphased_intervals
}


pub fn build_aware_contigs(target_sequences:&[Vec<u8>], target_intervals:&Vec<SeqInterval>, haplotypes:&FxHashMap<HaplotypeId,Haplotype>, min_length:usize) -> Vec<AwareContig> {
    
    let mut aware_contigs = retrieve_unphased_intervals(target_intervals, haplotypes, min_length)
        .into_iter()
        .map(|siv| AwareContig{
            interval: siv,
            is_separated: siv.beg != 0 || siv.end != target_sequences[siv.tid].len(),
            hid: None,
            sequence: target_sequences[siv.tid][siv.beg..siv.end].to_vec(),
            aligned_bases: 0
        }).collect_vec();

    aware_contigs.extend(haplotypes.values()
        .map(|ht| {
            let piv = ht.region();
            AwareContig {
                interval: piv,
                is_separated: true,
                hid: Some(ht.hid()),
                sequence: target_sequences[piv.tid][piv.beg..piv.end].to_vec(),
                aligned_bases: 0
            }
        })
    );

    aware_contigs.sort_unstable_by_key(|ctg| (ctg.interval(), ctg.hid()));
    aware_contigs
}


// aware_contigs should be sorted by interval
pub fn map_alignments_to_aware_contigs(seq_alignments: &[SeqAlignment], aware_contigs: &[AwareContig], seq2haplo: &FxHashMap<BamRecordId,Vec<HaplotypeId>>, ambiguous_reads:&FxHashSet<String>) -> Vec<AwareAlignment> {
    
    let mut aware_alignments: Vec<AwareAlignment> = vec![];

    for sa in seq_alignments.iter().filter(|sa| !ambiguous_reads.contains(sa.query_name())) {
        let sa_id = sa.record_id();

        let idx = aware_contigs.partition_point(|ctg| ctg.tid() < sa.tid() || (ctg.tid() == sa.tid() && ctg.end() <= sa.target_beg()));
        for (i,ctg) in aware_contigs[idx..].iter().take_while(|ctg| ctg.beg() < sa.target_end()).enumerate() {
            
            if ctg.is_phased() && (!seq2haplo.contains_key(&sa_id) || !seq2haplo[&sa_id].contains(&ctg.haplotype_id().unwrap())) {
                continue
            }

            let aware_id = idx+i;
            
            let aware_target_beg = std::cmp::max(sa.target_beg(), ctg.beg());
            let aware_target_end = std::cmp::min(sa.target_end(), ctg.end());
            
            let aware_ctg_beg = aware_target_beg - ctg.beg();
            let aware_ctg_end = aware_target_end - ctg.beg();
            let aware_ctg_range = (aware_ctg_beg, aware_ctg_end, ctg.length());

            let query_beg = if let Strand::Forward = sa.strand() { sa.query_beg() } else { sa.query_length() - sa.query_end() };
            let (mut aware_query_beg, aware_target_beg) = crate::seq::alignment::first_match_from(ctg.beg(), query_beg, sa.target_beg(), sa.cigar());
            let (mut aware_query_end, aware_target_end) = crate::seq::alignment::last_match_until(ctg.end(), query_beg, sa.target_beg(), sa.cigar());
            let aware_query_range = (aware_query_beg, aware_query_end, sa.query_length());

            if aware_query_end <= aware_query_beg || aware_target_end <= aware_target_beg {
                continue
            }

            let maptype = crate::seq::alignment::classify_mapping(aware_query_range, aware_ctg_range);

            if let Strand::Reverse = sa.strand() {
                (aware_query_beg, aware_query_end) = (sa.query_length() - aware_query_beg, sa.query_length() - aware_query_end);
            }

            if matches!(maptype, MappingType::QueryContained|MappingType::ReferenceContained|MappingType::DovetailPrefix|MappingType::DovetailSuffix) {
                aware_alignments.push(AwareAlignment{
                    query_name: sa.query_name().to_string(),
                    query_length: sa.query_length(),
                    query_beg: aware_query_beg,
                    query_end: aware_query_end,
                    strand: sa.strand(),
                    aware_id,
                    target_name: sa.target_name().to_string(),
                    target_beg: aware_target_beg,
                    target_end: aware_target_end,
                    mapq: sa.mapq(),
                    identity: sa.identity(),
                    query_aln_beg: sa.query_beg(),
                    query_aln_end: sa.query_end(),
                    mapping_type: maptype
                });
            }
        }
    }

    aware_alignments
}

/*
# TODO: use read_alignments instead of reference_alignments ?
# TODO: refactor code without the possibility of ambiguous mappings (break alignment paths, instead)
def map_sreads_to_aware_contigs(reference_alignments, aware_table, aware_intervals, segment_haplotypes, ambiguous_segments, allow_ambiguous=False):
    for reference_id, alignments in reference_alignments.items():
        for a in alignments:

            for aware_contig in (iv.data for iv in sorted(aware_intervals[reference_id].overlap(a.reference_start,a.reference_end))):
                
                if a.strand == '-': 
                    aware_qry_start, aware_qry_end = a.query_length-aware_qry_end, a.query_length-aware_qry_start

                if maptype in [MappingType.QUERY_CONTAINED, MappingType.REFERENCE_CONTAINED, MappingType.DOVETAIL_PREFIX, MappingType.DOVETAIL_SUFFIX]:
                    aware_alignment = AwareAlignment(a.query, a.query_length, aware_qry_start, aware_qry_end, a.strand, aware_contig.id, 
                                                     reference_id, aware_ref_start, aware_ref_end, a.mapq, a.identity(), a.query_start, a.query_end, maptype)
                    aware_alignments[a.query].append(aware_alignment)
    
    for query_id in ambiguous_reads:
        aware_alignments.pop(query_id,None)
    
    for query_id in aware_alignments:
        aware_alignments[query_id] = sorted(aware_alignments[query_id], key=lambda x:x.query_start)
        # if query_id in ['edge_8_517622-524797_h1_seq1']:
        #     logger.debug(f'---------- {query_id} ----------')
        #     for a in aware_alignments[query_id]:
        #         logger.debug(f'{a}')
        aware_alignments[query_id] = filter_alignments(aware_alignments[query_id],aware_table)
        # if query_id in ['edge_8_517622-524797_h1_seq1']:
        #     logger.debug(f'########## {query_id} ##########')
        #     for a in aware_alignments[query_id]:
        #         logger.debug(f'{a}')
    return aware_alignments
*/

/*
def filter_alignments(alignments, aware_table):
    sorted_aln = sorted(alignments, key=lambda x:(x.query_aln_start,x.query_start))
    filtered_aln = []
    for i,a in enumerate(sorted_aln):
        if a.query_end-a.query_start < 1:
            continue
        if a.maptype in [MappingType.DOVETAIL_PREFIX,MappingType.DOVETAIL_SUFFIX] and (0 < i < len(sorted_aln)-1):
            continue
        aware_contig = aware_table[a.aware_id]
        aware_aligned_bases = a.reference_end-a.reference_start
        aware_contig.alignedbases += aware_aligned_bases
        filtered_aln.append(a)
    return filtered_aln
*/

/*
def map_sreads_to_aware_contigs_old(reference_alignments, aware_table, aware_intervals, segment_haplotypes, ambiguous_segments, allow_ambiguous=False):
    aware_alignments = defaultdict(list)
    ambiguous_reads = set()
    for reference_id, alignments in reference_alignments.items():
        for a in alignments:
            segment_id = a.id()
            if (not allow_ambiguous) and (segment_id in ambiguous_segments):
                ambiguous_reads.add(a.query)
                continue
            for aware_contig in (iv.data for iv in sorted(aware_intervals[reference_id].overlap(a.reference_start,a.reference_end))):
                if aware_contig.phased and (aware_contig.name() not in segment_haplotypes[segment_id]): 
                    continue
                overlap_on_reference = overlap_length(a.reference_start, a.reference_end, aware_contig.reference_start, aware_contig.reference_end)
                aware_contig.alignedbases += overlap_on_reference
                query_start, reference_start = first_match_from(aware_contig.reference_start, a.query_start, a.reference_start, a.cigar)
                query_end, reference_end = last_match_until(aware_contig.reference_end, a.query_start, a.reference_start, a.cigar)
                if a.strand == '-': query_start,query_end = a.query_length-query_end, a.query_length-query_start
                query_aln_start, query_aln_end = (a.query_start,a.query_end) if a.strand == '+' else (a.query_length-a.query_end,a.query_length-a.query_start)
                aware_alignment = AwareAlignment(a.query, a.query_length, query_start, query_end, a.strand, aware_contig.id, reference_start, reference_end, a.mapq, a.identity(), query_aln_start, query_aln_end)
                aware_alignments[a.query].append((segment_id,aware_alignment))
    for query_id in ambiguous_reads:
        aware_alignments.pop(query_id,None)
    for query_id in aware_alignments:
        aware_alignments[query_id] = filter_alignments(aware_alignments[query_id],aware_table)
    return aware_alignments
*/