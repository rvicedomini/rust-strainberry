// use std::io::BufRead;
use std::sync::mpsc;
use std::thread;
use std::path::Path;

use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use anyhow::Result;
use itertools::Itertools;
use rust_htslib::bam::{self, Read,IndexedReader};
use tinyvec::ArrayVec;

use crate::cli::Options;
use crate::seq::SeqDatabase;


const BASES: [char;4] = ['A', 'C', 'G', 'T'];

#[derive(Debug, Clone)]
pub struct Var {
    pub tid: usize,
    pub pos: usize,
    pub depth: usize,
    pub alleles: ArrayVec<[(char,usize);2]>,
}

pub type VarDict = HashMap<usize,Vec<Var>>;
type VarPositions = HashMap<usize,HashSet<usize>>;


fn load_variants_at_positions_threaded(bam_path: &Path, ref_db: &SeqDatabase, positions: Option<&VarPositions>, opts: &Options) -> VarDict {

    let (tx, rx) = mpsc::channel();

    let index_to_tid = crate::bam::build_target_index(bam_path, ref_db);

    let ref_variants = if let Some(positions) = positions {
        positions.iter()
            .map(|(ref_idx,ref_positions)| (index_to_tid[*ref_idx], Some(ref_positions)))
            .collect_vec()
    } else {
        (0..ref_db.size())
            .map(|ref_idx| (index_to_tid[ref_idx], None))
            .collect_vec()
    };

    thread::scope(|scope| {
        for thread_id in 0..opts.nb_threads {
            let sender = tx.clone();
            let ref_variants_ref = &ref_variants;
            scope.spawn(move || {
                let mut bam_reader = IndexedReader::from_path(bam_path).unwrap();
                for i in (thread_id..ref_variants_ref.len()).step_by(opts.nb_threads) {
                    let (tid, target_positions) = ref_variants_ref[i];
                    let variants = load_variants_at_positions(&mut bam_reader, tid, target_positions, opts);
                    sender.send((tid,variants)).unwrap();
                }
            });
        }
    });

    (0..ref_variants.len())
        .flat_map(|_| rx.recv())
        .collect()
}


fn load_variants_at_positions(bam_reader: &mut IndexedReader, tid: usize, positions: Option<&HashSet<usize>>, opts: &Options) -> Vec<Var> {

    let mut variants: Vec<Var> = vec![];
    // let header = bam_reader.header().clone();

    bam_reader.fetch(tid as u32).unwrap();
    for pileup in bam_reader.pileup() {
        let pileup = pileup.unwrap();
        let pos = pileup.pos() as usize;

        if let Some(positions) = positions {
            if !positions.contains(&pos) {
                continue;
            }
        }

        let mut strand_counts = [0_usize; 8]; // A+,A-,C+,C-,G+,G-,T+,T-
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
                    let base = 0b11 & ((base >> 2) ^ (base >> 1));
                    let b = 2 * (base as usize) + (record.is_reverse() as usize);
                    strand_counts[b] += 1;
                }
            }
        }

        let depth = strand_counts.iter().sum::<usize>();
        let min_depth = std::cmp::max(opts.min_alt_count, (opts.min_alt_frac * (depth as f64)) as usize);
        let alleles: ArrayVec<[(usize,usize);4]> = strand_counts.chunks_exact(2)
            .map(|c| c[0]+c[1])
            .enumerate()
            .filter(|(_,cnt)| *cnt >= min_depth)
            .collect();

        if alleles.len() == 2 {

            // From longshot: https://github.com/pjedge/longshot/blob/99ace95bada7b360dc338deae65073590d6dc35d/src/main.rs#L478
            let b1 = alleles[0].0;
            let b2 = alleles[1].0;
            let counts = [
                strand_counts[2*b1] as u32, strand_counts[2*b1+1] as u32,
                strand_counts[2*b2] as u32, strand_counts[2*b2+1] as u32,
            ];

            let fishers_exact_pvalues = fishers_exact::fishers_exact(&counts)
                .expect("error calculating Fisher's exact test for strand bias.");

            if fishers_exact_pvalues.two_tail_pvalue >= opts.strand_bias_pvalue && counts.iter().all(|c| *c > 1) {
                // let tid = tid as usize;
                // let name = std::str::from_utf8(header.tid2name(tid as u32)).unwrap();
                // spdlog::trace!("{name} @ {pos}: {} | {:?}", fishers_exact_pvalues.two_tail_pvalue, counts);
                let alleles = alleles.into_iter().map(|(b,c)| (BASES[b], c)).collect();
                let var = Var{tid,pos,depth,alleles};
                variants.push(var);
            }

            // let alleles = alleles.into_iter().map(|(b,c)| (BASES[b], c)).collect();
            // let var = Var{tid,pos,depth,alleles};
            // variants.push(var);
        }
    }

    variants
}


pub fn load_variants_from_bam(bam_path:&Path, ref_db: &SeqDatabase, opts:&Options) -> VarDict {
    load_variants_at_positions_threaded(bam_path, ref_db, None, opts)
}


pub fn write_variants_to_file(path:&Path, variants:&VarDict, ref_db: &SeqDatabase) -> Result<()> {

    let mut writer = crate::utils::get_file_writer(path);
    // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
    for (ref_idx, ref_variants) in variants.iter() {
        for var in ref_variants {
            let ref_name = ref_db.names[*ref_idx].as_str();
            let var_pos = var.pos;
            let (var_nuc, _var_depth) = var.alleles[0];
            let (alt_nuc, _alt_depth) = var.alleles[1];
            let line = format!("{ref_name}\t{var_pos}\t{var_nuc}\t{alt_nuc}\n");
            writer.write_all(line.as_bytes())?;
        }
    }
    Ok(())
}


// pub fn load_variants_from_vcf(vcf_path:&Path, bam_path:&Path, ref_db: &SeqDatabase, opts:&Options) -> VarDict {

//     let vcf_reader = crate::utils::get_file_reader(vcf_path);
//     let positions = vcf_reader.lines()
//         .map_while(Result::ok)
//         .filter_map(|line| { // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
//             let line = line.trim();
//             if line.is_empty() || line.starts_with('#') { return None }
//             let mut record = line.split('\t');
//             let ref_name = record.next().unwrap();
//             let ref_index = ref_db.get_index(ref_name);
//             let pos = record.next().unwrap().parse::<usize>().unwrap();
//             let qual = record.nth(3).unwrap().parse::<usize>().unwrap();
//             if qual >= opts.min_var_qual { Some((ref_index,pos-1)) } else { None }
//         })
//         .sorted_unstable()
//         .collect_vec();

//     let positions = positions.into_iter()
//         .chunk_by(|&(ref_idx,_)| ref_idx).into_iter()
//         .map(|(ref_idx,ref_positions)| (
//             ref_idx,
//             ref_positions.into_iter().map(|(_,pos)| pos).collect()
//         ))
//         .collect::<VarPositions>();

//     load_variants_at_positions_threaded(bam_path, ref_db, Some(&positions), opts)
// }


// pub fn run_longcalld(reference_path: &Path, bam_path: &Path, ref_db: &SeqDatabase, vcf_path: &Path, opts: &Options) -> Result<VarDict> {

//     use std::process::{Command, Stdio};
//     use std::io::{BufRead,BufReader};

//     let lcd_preset = match opts.mode {
//         crate::cli::Mode::Hifi => "--hifi",
//         crate::cli::Mode::Nano => "--ont",
//     };
    
//     let args = ["call", lcd_preset, "-t", &opts.nb_threads.to_string(), reference_path.to_str().unwrap(), bam_path.to_str().unwrap()];
//     let stdout = Command::new("longcallD").args(args)
//         .stdout(Stdio::piped())
//         .stderr(Stdio::null())
//         .spawn().context("cannot run minimap2")?
//         .stdout
//         .with_context(|| "Could not capture standard output.")?;
    
//     let reader = BufReader::new(stdout);
//     let vcf_lines = reader.lines()
//         .map_while(Result::ok)
//         .filter_map(|line| {
//             let line = line.trim();
//             if !line.is_empty() { Some(line.to_string()) } else { None }
//         }).collect_vec();

//     let mut vcf_writer = crate::utils::get_file_writer(vcf_path);
//     for line in vcf_lines {
//         let fields = line.split('\t').collect_vec();
//         if line.starts_with('#') || (fields.len() >= 5 && fields[3].len() == 1 && fields[4].len() == 1) {
//             vcf_writer.write_all(line.as_bytes())?;
//             vcf_writer.write_all(b"\n")?;
//         }
//     }

//     let var_dict = load_variants_from_vcf(vcf_path, bam_path, ref_db, opts);
//     Ok(var_dict)
// }
