use std::io::BufRead;
use std::path::Path;

use rustc_hash::{FxHashMap,FxHashSet};

use rust_htslib::bam;
use rust_htslib::bam::Read;

use tinyvec::ArrayVec;

use crate::cli;
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

type VarPositions = FxHashSet<(usize,usize)>;

// load and filter variants (possibly restricting to a limited set of positions)
fn load_variants_at_positions(bam_path: &Path, positions: Option<&VarPositions>, opts: &cli::Options) -> VarDict {

    let mut bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header_view = bam_reader.header();

    let mut variants: VarDict = VarDict::default();
    
    for tid in 0..header_view.target_count() {
        bam_reader.fetch(tid).unwrap();

        for pileup in bam_reader.pileup() {
            let pileup = pileup.unwrap();
            let pos = pileup.pos() as usize;

            if let Some(positions) = positions {
                if !positions.contains(&(tid as usize,pos)) {
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
                    || record.is_supplementary()
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
                variants.entry(tid)
                    .or_insert(vec![])
                    .push(var);
            }
        }
    }

    variants.values_mut()
        .for_each(|vars| vars.sort_unstable_by_key(|v| v.pos));

    variants
}


pub fn load_variants_from_bam(bam_path:&Path, opts:&cli::Options) -> VarDict {
    load_variants_at_positions(bam_path, None, opts)
}


pub fn load_variants_from_vcf(vcf_path:&Path, bam_path:&Path, opts:&cli::Options) -> VarDict {

    let chrom2tid = utils::chrom2tid(bam_path);

    let vcf_reader = utils::get_file_reader(vcf_path);
    let positions: VarPositions = vcf_reader.lines().flatten()
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty() && !line.starts_with("#"))
        .map(|line| {
            let mut record = line.split('\t');
            let chrom = record.next().unwrap();
            let tid = *chrom2tid.get(chrom).unwrap();
            let pos = record.next().unwrap().parse::<usize>().unwrap();
            (tid,pos-1)
        })
        .collect();

    load_variants_at_positions(bam_path, Some(&positions), opts)
}



