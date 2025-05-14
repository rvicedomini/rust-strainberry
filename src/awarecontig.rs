use std::fmt;

use itertools::Itertools;
use ahash::AHashMap as HashMap;

use crate::phase::haplotype::{HaplotypeId, HaplotypeHit, Haplotype};
use crate::seq::SeqInterval;
use crate::alignment::{MappingType, SeqAlignment};
use crate::bam::BamRecordId;


#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    Haplotype(usize),
    #[default] Unphased,
    Read, // read strand
}

impl fmt::Display for SeqType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeqType::Haplotype(hid) => { write!(f, "Haplotype({hid})") }
            SeqType::Unphased => { write!(f, "Unphased") }
            SeqType::Read => { write!(f, "Read") }
        }
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct AwareContig {
    contig_type: SeqType,
    interval: SeqInterval,
    strand: u8,
    aligned_bases: usize,
}

impl AwareContig {

    pub fn new(contig_type:SeqType, interval:SeqInterval, strand:u8, aligned_bases:usize) -> Self {
        Self {
            contig_type,
            interval,
            strand,
            aligned_bases
        }
    }

    pub fn tid(&self) -> usize { self.interval.tid }
    pub fn beg(&self) -> usize { self.interval.beg }
    pub fn end(&self) -> usize { self.interval.end }
    pub fn length(&self) -> usize { self.interval.length() }
    pub fn strand(&self) -> u8 { self.strand }

    pub fn contig_type(&self) -> SeqType { self.contig_type }
    pub fn is_phased(&self) -> bool { matches!(self.contig_type, SeqType::Haplotype(_)) }
    
    pub fn hid(&self) -> Option<usize> {
        match self.contig_type {
            SeqType::Haplotype(hid) => Some(hid),
            _ => None
        }
    }

    pub fn haplotype_id(&self) -> Option<HaplotypeId> {
        Some(HaplotypeId{
            tid: self.interval.tid,
            beg: self.interval.beg,
            end: self.interval.end,
            hid: self.hid()?
        })
    }

    pub fn interval(&self) -> SeqInterval { self.interval }

    // pub fn phaseset_id(&self) -> SeqInterval { self.interval() }
    // pub fn region(&self) -> (usize,usize) { (self.beg(), self.end()) }

    pub fn depth(&self) -> f64 {
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
    pub aware_id: usize,
    pub query_idx: usize,
    pub query_len: usize,
    pub query_beg: usize,
    pub query_end: usize,
    pub strand: u8,
    pub target_idx: usize,
    pub target_beg: usize,
    pub target_end: usize,
    pub mapq: u8,
    pub identity: f64,
    pub query_aln_beg: usize,
    pub query_aln_end: usize,
    pub mapping_type: MappingType,
    pub dist: usize,
    pub nb_shared_snvs: usize,
    pub nb_alt_matches: usize,
}

impl AwareAlignment {

    pub fn is_ambiguous(&self) -> bool {
        self.nb_alt_matches > 0
    }

    pub fn target_interval(&self) -> SeqInterval {
        SeqInterval { tid: self.target_idx, beg: self.target_beg, end: self.target_end }
    }

}


fn retrieve_unphased_intervals(ref_intervals:&[SeqInterval], haplotypes:&HashMap<HaplotypeId,Haplotype>, min_length:usize) -> Vec<SeqInterval> {

    let mut unphased_intervals: Vec<SeqInterval> = vec![];
    let phased_intervals = haplotypes.values()
        .map(|ht| ht.region())
        .unique()
        .sorted_unstable()
        .collect_vec();

    for iv in ref_intervals {

        let first = phased_intervals.partition_point(|x| x.tid < iv.tid || (x.tid == iv.tid && x.beg < iv.beg));
        
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


pub fn build_aware_contigs(ref_intervals:&[SeqInterval], haplotypes:&HashMap<HaplotypeId,Haplotype>, min_length:usize) -> Vec<AwareContig> {
    
    let mut aware_contigs = retrieve_unphased_intervals(ref_intervals, haplotypes, min_length)
        .into_iter()
        .map(|interval| AwareContig{
            contig_type: SeqType::Unphased,
            interval,
            strand: b'+',
            aligned_bases: 0
        }).collect_vec();

    aware_contigs.extend(haplotypes.values()
        .map(|ht| {
            AwareContig {
                contig_type: SeqType::Haplotype(ht.hid()),
                interval: ht.region(),
                strand: b'+',
                aligned_bases: 0
            }
        })
    );

    aware_contigs.sort_by_key(|ctg| (ctg.interval(), ctg.hid()));
    spdlog::debug!("{} strain-aware contigs built", aware_contigs.len());

    aware_contigs
}


pub fn build_phased_contigs(haplotypes:&HashMap<HaplotypeId,Haplotype>) -> Vec<AwareContig> {
    
    let mut aware_contigs = haplotypes.values()
        .map(|ht| {
            AwareContig {
                contig_type: SeqType::Haplotype(ht.hid()),
                interval: ht.region(),
                strand: b'+',
                aligned_bases: 0
            }
        }).collect_vec();

    aware_contigs.sort_by_key(|ctg| (ctg.interval(), ctg.hid()));
    spdlog::debug!("{} strain-aware contigs built", aware_contigs.len());

    aware_contigs
}


pub fn map_sequences_to_aware_contigs(seq_alignments: &HashMap<usize,Vec<SeqAlignment>>, aware_contigs: &mut[AwareContig], seq2haplo: &HashMap<BamRecordId,Vec<HaplotypeHit>>) -> HashMap<usize,Vec<AwareAlignment>> {

    let mut read2aware = HashMap::new();
    for (read_idx, alignments) in seq_alignments {
        let aware_alignments = map_alignments_to_aware_contigs(alignments, aware_contigs, seq2haplo);
        let aware_alignments = filter_aware_alignments(aware_alignments, aware_contigs);
        read2aware.insert(*read_idx, aware_alignments);
    }

    read2aware
}


// aware_contigs should be sorted by interval
pub fn map_alignments_to_aware_contigs(alignments: &[SeqAlignment], aware_contigs: &[AwareContig], seq2haplo: &HashMap<BamRecordId,Vec<HaplotypeHit>>) -> Vec<AwareAlignment> {
    
    let mut aware_alignments: Vec<AwareAlignment> = vec![];

    for sa in alignments {
        let sa_id = sa.bam_record_id();

        let idx = aware_contigs.partition_point(|ctg| ctg.tid() < sa.tid() || (ctg.tid() == sa.tid() && ctg.end() <= sa.target_beg()));
        for (i,ctg) in aware_contigs[idx..].iter().take_while(|ctg| ctg.tid() == sa.tid() && ctg.beg() < sa.target_end()).enumerate() {

            if ctg.is_phased() && (!seq2haplo.contains_key(&sa_id) || !seq2haplo[&sa_id].iter().any(|hit| hit.hid == ctg.haplotype_id().unwrap())) {
                continue
            }

            let aware_id = idx+i;
            
            let aware_target_beg = std::cmp::max(sa.target_beg(), ctg.beg());
            let aware_target_end = std::cmp::min(sa.target_end(), ctg.end());
            
            let aware_ctg_beg = aware_target_beg - ctg.beg();
            let aware_ctg_end = aware_target_end - ctg.beg();
            let aware_ctg_range = (aware_ctg_beg, aware_ctg_end, ctg.length());

            let query_beg = if sa.is_forward() { sa.query_beg() } else { sa.query_length() - sa.query_end() };
            let (mut aware_query_beg, aware_target_beg) = crate::alignment::first_match_from(ctg.beg(), query_beg, sa.target_beg(), sa.cigar()).unwrap();
            let (mut aware_query_end, aware_target_end) = crate::alignment::last_match_until(ctg.end(), query_beg, sa.target_beg(), sa.cigar()).unwrap();
            let aware_query_range = (aware_query_beg, aware_query_end, sa.query_length());

            if aware_query_end <= aware_query_beg || aware_target_end <= aware_target_beg {
                continue
            }

            let maptype = crate::alignment::classify_mapping(aware_query_range, aware_ctg_range, 100, 0.5);

            if sa.is_reverse() {
                (aware_query_beg, aware_query_end) = (sa.query_length() - aware_query_end, sa.query_length() - aware_query_beg);
            }

            let mut nb_alt_matches = 0;
            let mut nb_shared_snvs = 0;
            let mut dist = 0;

            if ctg.is_phased() && seq2haplo.contains_key(&sa_id) {
                if let Some(hit) = seq2haplo[&sa_id].iter().find(|hit| hit.hid == ctg.haplotype_id().unwrap()) {
                    nb_alt_matches = hit.nb_alt;
                    nb_shared_snvs = hit.nb_pos;
                    dist = hit.dist;
                }
            }

            let awa = AwareAlignment{
                aware_id,
                query_idx: sa.query_index(),
                query_len: sa.query_length(),
                query_beg: aware_query_beg,
                query_end: aware_query_end,
                strand: sa.strand(),
                target_idx: sa.target_index(),
                target_beg: aware_target_beg,
                target_end: aware_target_end,
                mapq: sa.mapq(),
                identity: sa.identity(),
                query_aln_beg: sa.query_beg(),
                query_aln_end: sa.query_end(),
                mapping_type: maptype,
                dist,
                nb_shared_snvs,
                nb_alt_matches,
            };

            if (ctg.is_phased() && seq2haplo.contains_key(&sa_id)) || (!ctg.is_phased() && matches!(awa.mapping_type, MappingType::QueryContained|MappingType::ReferenceContained|MappingType::DovetailPrefix|MappingType::DovetailSuffix)) {
                aware_alignments.push(awa);
            }
        }
    }

    aware_alignments
}


pub fn filter_aware_alignments(mut aware_alignments: Vec<AwareAlignment>, aware_contigs: &mut [AwareContig]) -> Vec<AwareAlignment> {
    aware_alignments.sort_unstable_by_key(|a| (a.query_aln_beg,a.query_beg));
    aware_alignments.into_iter()
        .filter(|a| a.query_end - a.query_beg >= 1)
        .inspect(|a| {
            let ctg = &mut aware_contigs[a.aware_id];
            ctg.aligned_bases += a.target_end - a.target_beg;
        }).collect()
}