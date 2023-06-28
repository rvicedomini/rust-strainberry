use std::fs;
use std::path::Path;
use std::time::Instant;

use clap::Parser;

use rust_htslib::bam;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::Read;

use strainberry::opts::Options;
use strainberry::utils::*;

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

    for tid in 0..header_view.target_count() {
        bam_reader.fetch(tid).unwrap();
        let pileup = bam_reader.pileup();

        for p in pileup {
            let p = p.unwrap();

            let tid = pileup.tid() as usize;
            let chrom: String = target_names[tid].clone();
    }


    println!("Time: {}s | MaxRSS: {}GB", t_start.elapsed().as_secs_f64(), get_maxrss());
}
