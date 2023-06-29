use std::io::BufRead;
use std::path::Path;

use rustc_hash::{FxHashMap,FxHashSet};

use rust_htslib::bam;
use rust_htslib::bam::Read;

use tinyvec::ArrayVec;

use crate::opts::Options;
use crate::utils;


const BASES: [char;5] = ['A', 'C', 'G', 'T', 'N'];

#[derive(Debug, Clone)]
pub struct Var {
    pub tid: u32,
    pub pos: u32,
    pub depth: usize,
    pub alleles: ArrayVec<[(char,usize);5]>,
}

// load and filter variants at given positions
// if positions is empty, all pileups are scanned to find suitable variants
fn load_variants_at_positions(bam_path:&Path, positions:&FxHashSet<(u32,u32)>, opts:&Options) -> Vec<Var> {

    let mut bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header_view = bam_reader.header();

    let mut variants: Vec<Var> = Vec::new();
    
    for tid in 0..header_view.target_count() {
        bam_reader.fetch(tid).unwrap();

        for pileup in bam_reader.pileup() {
            let pileup = pileup.unwrap();
            let pos = pileup.pos();

            if !positions.is_empty() && !positions.contains(&(tid,pos)) {
                continue;
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

            if alleles.len() == 2 && alleles.iter().all(|&(nuc,_)| nuc != 'N') {
                variants.push(Var { tid, pos, depth, alleles });
            }
        }
    }

    variants
}


pub fn load_variants_from_bam(bam_path:&Path, opts:&Options) -> Vec<Var> {
    let positions = FxHashSet::default();
    load_variants_at_positions(bam_path, &positions, opts)
}


fn chrom2tid(bam_path:&Path) -> FxHashMap<String,u32> {
    let bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    bam_reader.header()
        .target_names().iter()
        .enumerate()
        .map(|(tid,&name)| (String::from_utf8_lossy(name).to_string(),tid as u32))
        .collect()
}


pub fn load_variants_from_vcf(vcf_path:&Path, bam_path:&Path, opts:&Options) -> Vec<Var> {

    let chrom2tid = chrom2tid(bam_path);

    let mut positions: FxHashSet<(u32, u32)> = FxHashSet::default();
    let vcf_reader = utils::get_file_reader(vcf_path);
    for line in vcf_reader.lines().flatten() {
        
        let line = line.trim();
        if line.is_empty() || line.starts_with("#") {
            continue
        }

        let mut record = line.split('\t');
        let chrom = record.next().unwrap();
        let tid = *chrom2tid.get(chrom).unwrap();
        let pos = record.next().unwrap().parse::<u32>().unwrap();
        assert!(pos > 0);
        positions.insert((tid,pos-1));
    }

    load_variants_at_positions(bam_path, &positions, opts)
}



