use std::io::BufRead;
use std::sync::mpsc;
use std::thread;
use std::path::Path;

use itertools::Itertools;

use rustc_hash::{FxHashMap,FxHashSet};

use rust_htslib::bam;
use rust_htslib::bam::{Read,IndexedReader};

use tinyvec::ArrayVec;

use crate::cli::Options;
use crate::utils;


const BASES: [char;5] = ['A', 'C', 'G', 'T', 'N'];

#[derive(Debug, Clone)]
pub struct Var {
    pub tid: usize,
    pub pos: usize,
    pub depth: usize,
    pub alleles: ArrayVec<[(char,usize);5]>,
}

pub type VarDict = FxHashMap<usize,Vec<Var>>;
type VarPositions = FxHashMap<usize,FxHashSet<usize>>;


fn load_variants_at_positions_threaded(bam_path: &Path, positions: Option<&VarPositions>, opts: &Options) -> VarDict {

    let (tx, rx) = mpsc::channel();

    let target_variants = if let Some(positions) = positions {
        positions.iter()
            .map(|(tid,target_positions)| (*tid,Some(target_positions)))
            .collect_vec()
    } else {
        utils::bam_target_intervals(bam_path)
            .into_iter()
            .map(|siv| (siv.tid, None))
            .collect_vec()
    };

    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let target_variants_ref = &target_variants;
            scope.spawn(move || {
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for i in (thread_id..target_variants_ref.len()).step_by(opts.nb_threads) {
                    let (tid, target_positions) = target_variants_ref[i];
                    let variants = load_variants_at_positions(&mut bam_reader, tid, target_positions, &opts);
                    sender.send((tid,variants)).unwrap();
                }
            });
        }
    });

    (0..target_variants.len())
        .flat_map(|_| rx.recv())
        .collect()
}


fn load_variants_at_positions(bam_reader: &mut IndexedReader, tid: usize, positions: Option<&FxHashSet<usize>>, opts: &Options) -> Vec<Var> {

    let mut variants: Vec<Var> = vec![];

    bam_reader.fetch(tid as u32).unwrap();
    for pileup in bam_reader.pileup() {
        let pileup = pileup.unwrap();
        let pos = pileup.pos() as usize;

        if let Some(positions) = positions {
            if !positions.contains(&pos) {
                continue;
            }
        }

        let mut counts = [0 as usize; 5]; // A,C,G,T,N
        for alignment in pileup.alignments() {

            let record = alignment.record();
            if record.mapq() < opts.min_mapq 
                || record.is_unmapped() 
                || record.is_secondary() 
                || record.is_quality_check_failed() 
                || record.is_duplicate()
            {
                continue;
            }

            if !alignment.is_del() && !alignment.is_refskip() {
                if let bam::pileup::Indel::None = alignment.indel() {
                    let base = alignment.record().seq()[alignment.qpos().unwrap()];
                    let b = match base {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => 4
                    };
                    counts[b] += 1;
                }
            }
        }
        
        let depth = counts.iter().sum::<usize>();
        let min_depth = std::cmp::max(opts.min_alt_count, (opts.min_alt_frac * (depth as f64)) as usize);
        let alleles: ArrayVec<[(char,usize);5]> = counts.iter()
            .enumerate()
            .filter(|&(_,cnt)| *cnt >= min_depth)
            .map(|(b,&cnt)| (BASES[b],cnt))
            .collect();

        if alleles.len() == 2 && alleles.iter().all(|&(base,_)| base != 'N') {
            let tid = tid as usize;
            let var = Var{tid,pos,depth,alleles};
            variants.push(var);
        }
    }

    variants
}


pub fn load_variants_from_bam(bam_path:&Path, opts:&Options) -> VarDict {
    load_variants_at_positions_threaded(bam_path, None, opts)
}


pub fn load_variants_from_vcf(vcf_path:&Path, bam_path:&Path, opts:&Options) -> VarDict {

    let bam_reader = bam::Reader::from_path(bam_path).unwrap();
    let bam_header = bam_reader.header();
    let chrom2tid = utils::chrom2tid(bam_header);

    let vcf_reader = utils::get_file_reader(vcf_path);
    let positions = vcf_reader.lines().flatten()
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty() && !line.starts_with("#"))
        .map(|line| {
            let mut record = line.split('\t');
            let chrom = record.next().unwrap();
            let tid = *chrom2tid.get(chrom).unwrap();
            let pos = record.next().unwrap().parse::<usize>().unwrap();
            (tid,pos-1)
        })
        .sorted_unstable()
        .collect_vec();

    let positions = positions.into_iter()
        .group_by(|&(tid,_)| tid)
        .into_iter()
        .map(|(tid,target_positions)| (
            tid,
            target_positions.into_iter().map(|(_,pos)| pos).collect()
        ))
        .collect::<VarPositions>();

    load_variants_at_positions_threaded(bam_path, Some(&positions), opts)
}



