use std::fs;
use std::path::Path;
use std::time::Instant;

use clap::Parser;

use strainberry::opts;
use strainberry::utils;
use strainberry::variant;


fn main() {
    
    let t_start = Instant::now();

    let opts = opts::Options::parse();

    let _fasta_path = Path::new(&opts.fasta_file);
    let bam_path = Path::new(&opts.bam_file);
    let output_dir = Path::new(&opts.output_dir);

    if fs::create_dir_all(output_dir).is_err() {
        println!("Cannot create output directory: \"{}\"", output_dir.display());
        std::process::exit(1);
    };

    let variants = if let Some(vcf_file) = &opts.vcf_file {
        variant::load_variants_from_vcf(Path::new(vcf_file), bam_path, &opts)
    } else {
        variant::load_variants_from_bam(bam_path, &opts)
    };
    
    for var in variants {
        println!("{var:?}");
    }

    println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
}
