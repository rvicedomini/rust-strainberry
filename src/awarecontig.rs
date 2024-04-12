use std::fmt;

use itertools::Itertools;
use rustc_hash::FxHashMap;

use crate::phase::haplotype::Haplotype;
use crate::seq::SeqInterval;
use crate::seq::alignment::MappingType;


#[derive(Debug, Clone)]
pub struct AwareContig {
    id: usize,
    interval: SeqInterval,
    is_separated: bool,
    haplotype_id: Option<usize>,
    sequence: Vec<u8>,
    aligned_bases: usize,
}

impl fmt::Display for AwareContig {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "AwareContig(id={}, iv={}, is_separated={}, ht_id={:?})", self.id, self.interval, self.is_separated, self.haplotype_id)
    }
}

impl AwareContig {

    fn is_phased(&self) -> bool {
        self.haplotype_id.is_some()
    }

    fn length(&self) -> usize {
        self.interval.end - self.interval.beg
    }

    fn interval(&self) -> SeqInterval {
        self.interval
    }

    fn phaseset_id(&self) -> SeqInterval {
        self.interval()
    }

    fn region(&self) -> (usize,usize) {
        (self.interval.beg, self.interval.end)
    }

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
struct AwareAlignment {
    query_id: String,
    query_length: usize,
    query_beg: usize,
    query_end: usize,
    strand: bool,
    aware_id: usize,
    reference_id: String,
    reference_beg: usize,
    reference_end: usize,
    mapq: u8,
    identity: f64,
    query_aln_beg: usize,
    query_aln_end: usize,
    mapping_type: MappingType
}


fn retrieve_unphased_intervals(target_intervals:&Vec<SeqInterval>, haplotypes:&FxHashMap<SeqInterval,Vec<Haplotype>>, min_length:usize) -> Vec<SeqInterval> {

    let mut unphased_intervals: Vec<SeqInterval> = vec![];
    let phased_intervals = haplotypes.keys().sorted_unstable().collect_vec();

    for iv in target_intervals {

        let first = phased_intervals.partition_point(|&x| x < iv);
        
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


pub fn build_aware_contigs(target_sequences:&[Vec<u8>], target_intervals:&Vec<SeqInterval>, haplotypes:&FxHashMap<SeqInterval,Vec<Haplotype>>, min_length:usize) -> Vec<AwareContig> {
    
    let mut aware_contigs = retrieve_unphased_intervals(target_intervals, haplotypes, min_length)
        .into_iter().enumerate()
        .map(|(aware_id,siv)| AwareContig{
            id: aware_id,
            interval: siv,
            is_separated: siv.beg != 0 || siv.end != target_sequences[siv.tid].len(),
            haplotype_id: None,
            sequence: target_sequences[siv.tid][siv.beg..siv.end].to_vec(),
            aligned_bases: 0
        }).collect_vec();

    aware_contigs.extend(haplotypes.values()
        .flatten()
        .zip(aware_contigs.len()..)
        .map(|(ht,aware_id)| {
            let piv = ht.region();
            AwareContig {
                id: aware_id,
                interval: piv,
                is_separated: true,
                haplotype_id: Some(ht.hid()),
                sequence: target_sequences[piv.tid][piv.beg..piv.end].to_vec(),
                aligned_bases: 0
            }
        })
    );

    aware_contigs
}



/*
# TODO: use read_alignments instead of reference_alignments ?
# TODO: refactor code without the possibility of ambiguous mappings (break alignment paths, instead)
def map_sreads_to_aware_contigs(reference_alignments, aware_table, aware_intervals, segment_haplotypes, ambiguous_segments, allow_ambiguous=False):
    aware_alignments = defaultdict(list)
    ambiguous_reads = set()
    for reference_id, alignments in reference_alignments.items():
        for a in alignments:
            segment_id = a.uid()
            if (not allow_ambiguous) and (segment_id in ambiguous_segments):
                ambiguous_reads.add(a.query)
                continue
            for aware_contig in (iv.data for iv in sorted(aware_intervals[reference_id].overlap(a.reference_start,a.reference_end))):
                
                if aware_contig.phased and all(ht_name != aware_contig.name() for ht_name,ht_dist in segment_haplotypes[segment_id]):
                    continue
                
                aware_ref_start = max(a.reference_start,aware_contig.reference_start)
                aware_ref_end = min(a.reference_end,aware_contig.reference_end)
                
                aware_ctg_start = aware_ref_start - aware_contig.reference_start
                aware_ctg_end = aware_ref_end - aware_contig.reference_start
                aware_ctg_range = (aware_ctg_start,aware_ctg_end,aware_contig.length())

                query_start = a.query_start if a.strand == '+' else (a.query_length-a.query_end)
                aware_qry_start, aware_ref_start = first_match_from(aware_contig.reference_start, query_start, a.reference_start, a.cigar)
                aware_qry_end, aware_ref_end = last_match_until(aware_contig.reference_end, query_start, a.reference_start, a.cigar)
                aware_qry_range = (aware_qry_start,aware_qry_end,a.query_length)

                if aware_qry_end-aware_qry_start <= 0 or aware_ref_end-aware_ref_start <= 0:
                    continue

                maptype = classify_mapping(aware_qry_range, aware_ctg_range)

                if a.strand == '-': 
                    aware_qry_start, aware_qry_end = a.query_length-aware_qry_end, a.query_length-aware_qry_start

                if maptype in [MappingType.QUERY_CONTAINED, MappingType.REFERENCE_CONTAINED, MappingType.DOVETAIL_PREFIX, MappingType.DOVETAIL_SUFFIX]:
                    aware_alignment = AwareAlignment(a.query, a.query_length, aware_qry_start, aware_qry_end, a.strand, aware_contig.id, 
                                                     reference_id, aware_ref_start, aware_ref_end, a.mapq, a.identity(), a.query_start, a.query_end, maptype)
                    aware_alignments[a.query].append(aware_alignment)

                # if a.query in ['edge_4_10631-10663_h1_seq1']:
                #     logger.debug(f'---------- {segment_id} ----------')
                #     logger.debug(f'aware_ctg = {aware_contig.name()}')
                #     logger.debug(f'aware_qry_range = {aware_qry_range}')
                #     logger.debug(f'aware_ref_interval = ({aware_ref_start},{aware_ref_end})')
                #     logger.debug(f'aware_range = {aware_ctg_range}')
                #     logger.debug(f'{str(maptype)}')
    
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