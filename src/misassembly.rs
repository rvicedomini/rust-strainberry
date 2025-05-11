use std::sync::mpsc;
use std::thread;
use std::path::Path;

use itertools::Itertools;
use rust_htslib::bam::{Reader, Read, IndexedReader, HeaderView};
use rust_htslib::bam::record::Aux;

use crate::cli::Options;
use crate::seq::{SeqDatabase, SeqInterval};
use crate::alignment::SeqAlignment;


#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Misassembly(usize,usize); // target_index, target_position


#[derive(Debug)]
struct AlignedBlock {
    query_beg: usize,
    query_end: usize,
    strand: u8,
    target_idx: usize,
    target_beg: usize,
    target_end: usize,
}


pub fn partition_reference(bam_path: &Path, ref_db: &SeqDatabase, read_db: &SeqDatabase, opts: &Options) -> Vec<SeqInterval> {

    let (tx, rx) = mpsc::channel();

    let bam_intervals = {
        let mut bam_intervals = crate::bam::bam_intervals(bam_path);
        bam_intervals.sort_unstable_by(|a,b| a.length().cmp(&b.length()).reverse());
        bam_intervals
    };

    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let bam_intervals_ref = &bam_intervals;
            scope.spawn(move || {
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for i in (thread_id..bam_intervals_ref.len()).step_by(opts.nb_threads) {
                    let candidates = find_misassemblies(&mut bam_reader, &bam_intervals_ref[i], ref_db, read_db, opts);
                    sender.send(candidates).unwrap();
                }
            });
        }
    });

    let candidates = (0..bam_intervals.len())
        .flat_map(|_| rx.recv().unwrap())
        .collect_vec();

    let bam_reader = Reader::from_path(bam_path).unwrap();
    let bam_header = bam_reader.header();
    cluster_misassemblies(candidates, bam_header, opts)
    
}


fn find_misassemblies(bam_reader: &mut IndexedReader, region: &SeqInterval, ref_db: &SeqDatabase, read_db: &SeqDatabase, opts: &Options) -> Vec<Misassembly> {

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

        let seqalign = SeqAlignment::from_bam_record(&record, &bam_header, ref_db, read_db);

        for (a,b) in seqalign.aligned_blocks().tuple_windows() {
            let [_a_query_beg,a_query_end,_a_target_beg,a_target_end] = a;
            let [b_query_beg,_b_query_end,b_target_beg,_b_target_end] = b;
            if b_target_beg - a_target_end >= opts.min_indel {
                candidates.push(Misassembly(seqalign.target_index(), a_target_end));
                candidates.push(Misassembly(seqalign.target_index(), b_target_beg));
            } else if b_query_beg - a_query_end >= opts.min_indel {
                candidates.push(Misassembly(seqalign.target_index(), b_target_beg));
            }
        }

        if !record.is_supplementary() {

            let mut read_alignments = vec![AlignedBlock{
                query_beg: seqalign.query_beg(),
                query_end: seqalign.query_end(),
                strand: seqalign.strand(),
                target_idx: seqalign.target_index(),
                target_beg: seqalign.target_beg(),
                target_end: seqalign.target_end(),
            }];

            if let Ok(Aux::String(supplementary_alignments)) = record.aux(b"SA") {

                let supplementary_alignments = supplementary_alignments
                    .split(';')
                    .take_while(|sa| !sa.is_empty())
                    .map(|sa| sa.split(',').collect_vec())
                    .filter(|sa| sa[4].parse::<u8>().unwrap() >= opts.min_mapq)
                    .collect_vec();

                for sa in supplementary_alignments {
                    let target_name = sa[0];
                    let target_idx = ref_db.get_index(target_name);
                    let target_beg = sa[1].parse::<usize>().unwrap() - 1;
                    let strand = sa[2].as_bytes()[0];
                    assert!(strand == b'+' || strand == b'-');
                    let cigar = crate::bam::parse_cigar_bytes(sa[3].as_bytes());
                    let [mut query_beg,mut query_end,target_beg,target_end] = crate::bam::intervals_from_cigar(&cigar, target_beg, 0);
                    let query_length = crate::bam::seq_length_from_cigar(&cigar, true);
                    if strand == b'-' {
                        (query_beg,query_end) = (query_length-query_end, query_length-query_beg);
                    }
                    read_alignments.push(AlignedBlock{ query_beg, query_end, strand, target_idx, target_beg, target_end });
                }
            }

            read_alignments.sort_unstable_by_key(|t| t.query_beg);
            candidates.append(&mut misassemblies_from_alignments(&read_alignments, seqalign.query_length(), opts));
        }
    }

    candidates
}


fn misassemblies_from_alignments(alignments: &[AlignedBlock], read_length: usize, opts: &Options) -> Vec<Misassembly> {
    
    let mut candidates = vec![];
    if alignments.is_empty() {
        return candidates;
    }

    let first = alignments.first().unwrap();
    if first.query_beg >= opts.min_overhang {
        candidates.push(if first.strand == b'+' {
            Misassembly(first.target_idx, first.target_beg)
        } else {
            Misassembly(first.target_idx, first.target_end)
        });
    }

    for w in alignments.windows(2) {
        let (a,b) = (&w[0],&w[1]);
        candidates.push(if a.strand == b'+' {
            Misassembly(a.target_idx, a.target_end)
        } else {
            Misassembly(a.target_idx, a.target_beg)
        });
        candidates.push(if b.strand == b'+' {
            Misassembly(b.target_idx, b.target_beg)
        } else {
            Misassembly(b.target_idx, b.target_end)
        });
    }

    let last = alignments.last().unwrap();
    if read_length - last.query_end >= opts.min_overhang {
        candidates.push(if last.strand == b'+' {
            Misassembly(last.target_idx, last.target_end)
        } else {
            Misassembly(last.target_idx, last.target_beg)
        });
    }
    
    candidates
}


fn cluster_misassemblies(mut candidates: Vec<Misassembly>, bam_header: &HeaderView, opts: &Options) -> Vec<SeqInterval> {

    candidates.sort_unstable();

    let mut clusters = vec![];
    for ma in candidates {
        if clusters.is_empty() {
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

    let mut intervals = vec![];
    for (tid, target_clusters) in &clusters.into_iter().filter(|clust| clust.len() >= opts.min_alt_count).chunk_by(|clust| clust[0].0) {
        let target_length = bam_header.target_len(tid as u32).unwrap() as usize;
        let mut target_pos = 0;
        for clust in target_clusters {
            let Misassembly(_, me_beg) = clust.first().unwrap();
            let Misassembly(_, me_end) = clust.last().unwrap();
            let Misassembly(_, me_median_pos) = &clust[clust.len()/2];
            // TODO: try splitting simply at median position, removing target[beg:end] could be too much?
            if *me_median_pos >= opts.min_overhang && target_length-me_median_pos >= opts.min_overhang {
                intervals.push(SeqInterval{ tid, beg: target_pos, end: *me_beg });
                target_pos = *me_end;
            }
        }
        intervals.push(SeqInterval { tid, beg: target_pos, end: target_length });
    }

    intervals
}
