use std::fmt;

use itertools::Itertools;
use rustc_hash::{FxHashMap,FxHashSet};

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


pub fn map_succinct_seqs_to_aware_contigs(seq_alignments: &FxHashMap<String,Vec<SeqAlignment>>, aware_contigs: &mut[AwareContig], seq2haplo: &FxHashMap<BamRecordId,Vec<HaplotypeId>>, ambiguous_reads:&FxHashSet<String>) -> FxHashMap<String,Vec<AwareAlignment>> {

    let mut read2aware: FxHashMap<String,Vec<AwareAlignment>> = FxHashMap::default();
    for (name, alignments) in seq_alignments {
        if ambiguous_reads.contains(name) {
            continue
        }
        let aware_alignments = map_alignments_to_aware_contigs(alignments, aware_contigs, seq2haplo);
        let aware_alignments = filter_aware_alignments(aware_alignments, aware_contigs);
        read2aware.insert(name.clone(), aware_alignments);
    }

    read2aware
}


// aware_contigs should be sorted by interval
pub fn map_alignments_to_aware_contigs(seq_alignments: &[SeqAlignment], aware_contigs: &[AwareContig], seq2haplo: &FxHashMap<BamRecordId,Vec<HaplotypeId>>) -> Vec<AwareAlignment> {
    
    let mut aware_alignments: Vec<AwareAlignment> = vec![];

    for sa in seq_alignments {
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
            let (mut aware_query_beg, aware_target_beg) = crate::seq::alignment::first_match_from(ctg.beg(), query_beg, sa.target_beg(), sa.cigar()).unwrap();
            let (mut aware_query_end, aware_target_end) = crate::seq::alignment::last_match_until(ctg.end(), query_beg, sa.target_beg(), sa.cigar()).unwrap();
            let aware_query_range = (aware_query_beg, aware_query_end, sa.query_length());

            if aware_query_end <= aware_query_beg || aware_target_end <= aware_target_beg {
                continue
            }

            let maptype = crate::seq::alignment::classify_mapping(aware_query_range, aware_ctg_range, 100, 0.05);

            if let Strand::Reverse = sa.strand() {
                (aware_query_beg, aware_query_end) = (sa.query_length() - aware_query_end, sa.query_length() - aware_query_beg);
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


pub fn filter_aware_alignments(mut aware_alignments: Vec<AwareAlignment>, aware_contigs: &mut [AwareContig]) -> Vec<AwareAlignment> {
    let nb_aware_alignments = aware_alignments.len();
    let mut filtered_alignments = Vec::with_capacity(nb_aware_alignments);
    
    aware_alignments.sort_unstable_by_key(|a| (a.query_aln_beg,a.query_beg));
    for (i,a) in aware_alignments.into_iter().enumerate() {
        if a.query_end - a.query_beg < 1 {
            continue
        }
        if matches!(a.mapping_type, MappingType::DovetailPrefix|MappingType::DovetailSuffix) && 0 < i && i < nb_aware_alignments-1 {
            continue
        }
        let aware_contig = &mut aware_contigs[a.aware_id];
        aware_contig.aligned_bases += a.target_end - a.target_beg;
        filtered_alignments.push(a);
    }

    filtered_alignments
}