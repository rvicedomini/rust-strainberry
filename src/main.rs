use std::fs;
use std::path::Path;
use std::time::Instant;

use clap::Parser;

use rust_htslib::bam;
use rust_htslib::bam::Read;
// use rust_htslib::bam::pileup::Indel;

use itertools::Itertools;

use tinyvec::ArrayVec;

use strainberry::opts::Options;
use strainberry::utils::*;


const BASES: [char;5] = ['A', 'C', 'G', 'T', 'N'];


fn main() {
    
    let t_start = Instant::now();

    let opts = Options::parse();

    let fasta_path = Path::new(&opts.fasta);
    let bam_path = Path::new(&opts.bam);
    let output_dir = Path::new(&opts.output_dir);

    if fs::create_dir_all(output_dir).is_err() {
        println!("Cannot create output directory: \"{}\"", output_dir.display());
        std::process::exit(1);
    };

    let mut bam_reader = bam::IndexedReader::from_path(bam_path).unwrap();
    let header_view = bam_reader.header();
    let target_names = header_view.target_names().iter().map(|s| String::from_utf8_lossy(s).to_string()).collect_vec();

    for tid in 0..header_view.target_count() {
        bam_reader.fetch(tid).unwrap();
        println!("Processing: {}", target_names[tid as usize]);

        for pileup in bam_reader.pileup() {
            let pileup = pileup.unwrap();
            // println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());

            let mut depth: usize = 0;
            let mut counts = [0 as usize; 5]; // A,C,G,T,N

            // TODO: consider a short gap as an unknown character
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
                        depth += 1;
                    }
                }
            }

            let mut indices = [0,1,2,3,4];
            indices.sort_unstable_by_key(|i| std::cmp::Reverse(counts[*i]));

            // println!("Counts: {counts:?}, Indices: {indices:?}");

            let dom_base = BASES[indices[0]];
            let dom_count = counts[indices[0]];

            let alt_base = BASES[indices[1]];
            let alt_count = counts[indices[1]];

            let dom_frac: f64 = (dom_count as f64) / (depth as f64);
            let alt_frac: f64 = (alt_count as f64) / (depth as f64);

            if alt_base == 'N' || alt_count < opts.min_alt_count || alt_frac < opts.min_alt_frac {
                continue;
            }

            println!("Candidate SNV at {}:{} {dom_base}({dom_count}/{dom_frac:.2}) | {alt_base}({alt_count}/{alt_frac:.2})", target_names[tid as usize], pileup.pos());
        }
    }


    println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), get_maxrss());
}
