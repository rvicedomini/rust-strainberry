use std::io::BufRead;
use std::sync::mpsc;
use std::thread;
use std::path::Path;

use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use itertools::Itertools;
use rust_htslib::bam::{self, Read,IndexedReader};
use tinyvec::ArrayVec;

use crate::cli::Options;


const BASES: [char;5] = ['A', 'C', 'G', 'T', 'N'];

#[derive(Debug, Clone)]
pub struct Var {
    pub tid: usize,
    pub pos: usize,
    pub depth: usize,
    pub alleles: ArrayVec<[(char,usize);5]>,
}

pub type VarDict = HashMap<usize,Vec<Var>>;
type VarPositions = HashMap<usize,HashSet<usize>>;


fn load_variants_at_positions_threaded(bam_path: &Path, positions: Option<&VarPositions>, opts: &Options) -> VarDict {

    let (tx, rx) = mpsc::channel();

    let target_variants = if let Some(positions) = positions {
        positions.iter()
            .map(|(tid,target_positions)| (*tid,Some(target_positions)))
            .collect_vec()
    } else {
        crate::bam::bam_target_intervals(bam_path)
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
                    let variants = load_variants_at_positions(&mut bam_reader, tid, target_positions, opts);
                    sender.send((tid,variants)).unwrap();
                }
            });
        }
    });

    (0..target_variants.len())
        .flat_map(|_| rx.recv())
        .collect()
}


fn load_variants_at_positions(bam_reader: &mut IndexedReader, tid: usize, positions: Option<&HashSet<usize>>, opts: &Options) -> Vec<Var> {

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

        let mut counts = [0_usize; 5]; // A,C,G,T,N
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
            // let tid = tid as usize;
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
    let chrom2tid = crate::bam::chrom2tid(bam_header);

    let vcf_reader = crate::utils::get_file_reader(vcf_path);
    let positions = vcf_reader.lines()
        .map_while(Result::ok)
        .filter_map(|line| { // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') { return None }
            let mut record = line.split('\t');
            let chrom = record.next().unwrap();
            let tid = *chrom2tid.get(chrom).unwrap();
            let pos = record.next().unwrap().parse::<usize>().unwrap();
            let qual = record.nth(3).unwrap().parse::<usize>().unwrap();
            if qual >= opts.min_var_qual { Some((tid,pos-1)) } else { None }
        })
        .sorted_unstable()
        .collect_vec();

    let positions = positions.into_iter()
        .chunk_by(|&(tid,_)| tid)
        .into_iter()
        .map(|(tid,target_positions)| (
            tid,
            target_positions.into_iter().map(|(_,pos)| pos).collect()
        ))
        .collect::<VarPositions>();

    load_variants_at_positions_threaded(bam_path, Some(&positions), opts)
}


// From longshot: https://github.com/pjedge/longshot/blob/99ace95bada7b360dc338deae65073590d6dc35d/src/main.rs#L478
// use Fishers exact test to filter variants for which allele observations are biased toward one strand or the other
pub fn filter_variants(vardict:HashMap<usize,Vec<Var>>) -> HashMap<usize,Vec<Var>> {

    let mut _retained = HashMap::new();

    for var in vardict.into_values().flatten() {
        if var.alleles.len() != 2 {
            continue
        }

        // let counts: [u32; 4] = [
        //     var.allele_counts_forward[0] as u32,
        //     var.allele_counts_reverse[0] as u32,
        //     var.allele_counts_forward[1] as u32,
        //     var.allele_counts_reverse[1] as u32,
        // ];
        // let fishers_exact_pvalues = fishers_exact(&counts)
        //     .chain_err(|| "Error calculating Fisher's exact test for strand bias.")?;

        // //println!("{:?} {:?} {:?}  {:?}",&counts, fishers_exact_pvalues.two_tail_pvalue, fishers_exact_pvalues.less_pvalue, fishers_exact_pvalues.greater_pvalue);
        // var.strand_bias_pvalue = if fishers_exact_pvalues.two_tail_pvalue <= 500.0 {
        //     *PHREDProb::from(Prob(fishers_exact_pvalues.two_tail_pvalue))
        // } else {
        //     500.0
        // };

        // if fishers_exact_pvalues.two_tail_pvalue < strand_bias_pvalue_cutoff {
        //     var.filter.add_filter(VarFilter::StrandBias);
        //     var.genotype = Genotype(0, 0);
        //     var.gq = 0.0;
        // }
    }

    // for f in 0..flist.len() {
    //     flist[f].calls.retain(|&c| {
    //         !varlist.lst[c.var_ix as usize]
    //             .filter
    //             .has_filter(VarFilter::StrandBias)
    //     });
    // }

    _retained
}
