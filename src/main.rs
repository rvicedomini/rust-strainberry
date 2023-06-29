use std::fs;
use std::path::Path;
use std::time::Instant;

use clap::Parser;

use rustc_hash::FxHashMap;
use strainberry::cli;
use strainberry::utils;
use strainberry::variant;

use needletail::{parse_fastx_file, Sequence};


fn main() {
    
    let t_start = Instant::now();

    let opts = cli::Options::parse();

    let fasta_path = Path::new(&opts.fasta_file);
    let bam_path = Path::new(&opts.bam_file);
    let output_dir = Path::new(&opts.output_dir);

    if fs::create_dir_all(output_dir).is_err() {
        println!("Cannot create output directory: \"{}\"", output_dir.display());
        std::process::exit(1);
    };

    let variants = if let Some(vcf_file) = &opts.vcf_file {
        println!("Loading and filtering variants from: {vcf_file}");
        variant::load_variants_from_vcf(Path::new(vcf_file), bam_path, &opts)
    } else {
        println!("Looking for potential variants");
        variant::load_variants_from_bam(bam_path, &opts)
    };

    println!("Found {} variants", variants.len());

    println!("Loading sequences from: {}", fasta_path.canonicalize().unwrap().display());

    let mut reference_sequences: FxHashMap<String,Vec<u8>> = FxHashMap::default();
    let mut fasta_reader = parse_fastx_file(fasta_path).expect("Cannot open reference file");
    while let Some(record) = fasta_reader.next() {
        let record = record.unwrap();
        let sid = String::from_utf8_lossy(record.id()).to_string();
        let sequence = record.normalize(false).to_vec();
        reference_sequences.insert(sid, sequence);
    }

    println!("Loaded {} sequences", reference_sequences.len());

    println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
}
