use std::sync::mpsc;
use std::thread;
use std::path::Path;

use itertools::Itertools;
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::{Reader, Read, IndexedReader};
use rust_htslib::bam::record::Aux;

use crate::cli::Options;
use crate::seq::SeqInterval;
use crate::seq::alignment::{Strand,SeqAlignment};
use crate::utils;


#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Misassembly(String,usize);


#[derive(Debug)]
struct AlignedBlock {
    query_beg: usize,
    query_end: usize,
    strand: Strand,
    target_name: String,
    target_beg: usize,
    target_end: usize,
}


pub fn partition_reference(bam_path: &Path, opts: &Options) -> Vec<SeqInterval> {

    let (tx, rx) = mpsc::channel();
    
    let target_intervals = utils::bam_target_intervals(bam_path)
        .into_iter()
        .sorted_unstable_by_key(|iv| iv.end-iv.beg)
        .collect_vec();

    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let target_intervals_ref = &target_intervals;
            scope.spawn(move || {
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for i in (thread_id..target_intervals_ref.len()).step_by(opts.nb_threads) {
                    let intervals = find_misassemblies(&mut bam_reader, &target_intervals_ref[i], &opts);
                    sender.send(intervals).unwrap();
                }
            });
        }
    });

    let candidates = (0..target_intervals.len())
        .flat_map(|_| rx.recv().unwrap())
        .collect_vec();

    let bam_reader = Reader::from_path(bam_path).unwrap();
    let bam_header = bam_reader.header();
    cluster_misassemblies(candidates, &bam_header, opts)
    
}


fn find_misassemblies(bam_reader: &mut IndexedReader, region: &SeqInterval, opts: &Options) -> Vec<Misassembly> {

    let mut candidates = vec![];
    let bam_header = bam_reader.header().clone();

    let fetch_definition = (region.tid as u32, region.beg as i64, region.end as i64);
    bam_reader.fetch(fetch_definition).expect("Failed fetching records from BAM");
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

        let aligned_blocks = seqalign.aligned_blocks().collect_vec();
        for w in aligned_blocks.windows(2) {
            let [_a_query_beg,a_query_end,_a_target_beg,a_target_end] = w[0];
            let [b_query_beg,_b_query_end,b_target_beg,_b_target_end] = w[1];
            if b_target_beg - a_target_end >= opts.min_indel {
                candidates.push(Misassembly(seqalign.target_name().to_string(),a_target_end));
                candidates.push(Misassembly(seqalign.target_name().to_string(),b_target_beg));
            } else if b_query_beg - a_query_end >= opts.min_indel {
                candidates.push(Misassembly(seqalign.target_name().to_string(),b_target_beg));
            }
        }

        if !record.is_supplementary() {

            let mut read_alignments = vec![AlignedBlock{
                query_beg: seqalign.query_beg(),
                query_end: seqalign.query_end(),
                strand: seqalign.strand(),
                target_name: seqalign.target_name().to_string(), 
                target_beg: seqalign.target_beg(),
                target_end: seqalign.target_end(),
            }];

            if let Some(Aux::String(supplementary_alignments)) = record.aux(b"SA").ok() {

                let supplementary_alignments = supplementary_alignments
                    .split(";")
                    .take_while(|sa| sa.len() > 0)
                    .map(|sa| sa.split(",").collect_vec())
                    .filter(|sa| sa[4].parse::<u8>().unwrap() >= opts.min_mapq)
                    .collect_vec();

                for sa in supplementary_alignments {
                    let target_name = sa[0].to_string();
                    let target_pos = sa[1].parse::<usize>().unwrap() - 1;
                    let strand = Strand::from_str(sa[2]).unwrap();
                    let cigar = utils::parse_cigar_bytes(sa[3].as_bytes());
                    let [mut query_beg,mut query_end,target_beg,target_end] = utils::intervals_from_cigar(&cigar, target_pos);
                    let query_length = utils::seq_length_from_cigar(&cigar, true);
                    if matches!(strand, Strand::Reverse) {
                        (query_beg,query_end) = (query_length-query_end, query_length-query_beg);
                    }
                    read_alignments.push(AlignedBlock{ query_beg, query_end, strand, target_name, target_beg, target_end });
                }
            }

            read_alignments.sort_unstable_by_key(|t| t.query_beg);
            candidates.append(&mut misassemblies_from_alignments(&read_alignments, seqalign.query_length(), opts));
        }
    }

    candidates
}


fn misassemblies_from_alignments(alignments: &Vec<AlignedBlock>, read_length: usize, opts: &Options) -> Vec<Misassembly> {
    
    let mut candidates = vec![];
    if alignments.len() == 0 {
        return candidates;
    }

    let first = alignments.first().unwrap();
    if first.query_beg >= opts.min_overhang {
        candidates.push(if let Strand::Forward = first.strand {
            Misassembly(first.target_name.clone(), first.target_beg)
        } else {
            Misassembly(first.target_name.clone(), first.target_end)
        });
    }

    for w in alignments.windows(2) {
        let (a,b) = (&w[0],&w[1]);
        candidates.push(if let Strand::Forward = a.strand {
            Misassembly(a.target_name.clone(), a.target_end)
        } else {
            Misassembly(a.target_name.clone(), a.target_beg)
        });
        candidates.push(if let Strand::Forward = b.strand {
            Misassembly(b.target_name.clone(), b.target_beg)
        } else {
            Misassembly(b.target_name.clone(), b.target_end)
        });
    }

    let last = alignments.last().unwrap();
    if read_length - last.query_end >= opts.min_overhang {
        candidates.push(if let Strand::Forward = last.strand {
            Misassembly(last.target_name.clone(), last.target_end)
        } else {
            Misassembly(last.target_name.clone(), last.target_beg)
        });
    }
    
    candidates
}


fn cluster_misassemblies(mut candidates: Vec<Misassembly>, bam_header: &HeaderView, opts: &Options) -> Vec<SeqInterval> {

    candidates.sort_unstable();

    let mut clusters = vec![];
    for ma in candidates {
        if clusters.len() == 0 {
            clusters.push(vec![ma]);
            continue
        }

        let cluster = clusters.last_mut().unwrap();
        
        if cluster[0].0 != ma.0 || cluster.iter().all(|x| ma.1 - x.1 > 30) { // TODO: get the "30" as parameter from opts
            clusters.push(vec![ma]);
        } else {
            cluster.push(ma);
        }
    }

    
    let chrom2tid = utils::chrom2tid(bam_header);

    let mut intervals = vec![];
    for (tid, target_clusters) in &clusters.into_iter().filter(|clust| clust.len() >= opts.min_alt_count).group_by(|clust| chrom2tid[&clust[0].0]) {
        let target_length = bam_header.target_len(tid as u32).unwrap() as usize;
        let mut target_pos = 0;
        for clust in target_clusters {
            let Misassembly(_, me_beg) = clust.first().unwrap();
            let Misassembly(_, me_end) = clust.last().unwrap();
            let Misassembly(_, me_median_pos) = &clust[clust.len()/2];
            if *me_median_pos >= opts.min_overhang && target_length-me_median_pos >= opts.min_overhang {
                intervals.push(SeqInterval{ tid, beg: target_pos, end: *me_beg });
                target_pos = *me_end;
            }
        }
        intervals.push(SeqInterval { tid, beg: target_pos, end: target_length });
    }

    intervals
}
