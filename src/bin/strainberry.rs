use std::fs;
use std::path::Path;
use std::time::Instant;

use clap::Parser;
use rustc_hash::FxHashMap;

use strainberry::cli;
use strainberry::misassembly;
use strainberry::phase;
use strainberry::seq::read::load_bam_sequences;
use strainberry::utils;
use strainberry::variant;


fn main() {
    
    let t_start = Instant::now();

    let opts = cli::Options::parse();

    let fasta_path = Path::new(&opts.fasta_file);
    let bam_path = Path::new(&opts.bam_file);
    let output_dir = Path::new(&opts.output_dir);

    // TODO?
    // create output directory in cli module ? (after option validation)
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
    println!("  {} variants found", variants.values().map(|vars| vars.len()).sum::<usize>());

    println!("Loading sequences from: {}", fasta_path.canonicalize().unwrap().display());
    let target_names = utils::bam_target_names(bam_path);
    let target_sequences = utils::load_sequences(fasta_path, bam_path);
    println!("  {} sequences loaded", target_sequences.len());
    
    // TODO:
    // Consider estimating lookback length when a flag "--auto-lookback" is provided
    // utils::estimate_lookback(bam_path, 1000)
    println!("Lookback {} bp", opts.lookback);

    println!("Splitting reference at putative misjoins");
    let target_intervals = misassembly::partition_reference(bam_path, &opts);

    println!("Loading reads from BAM");
    let read_sequences: FxHashMap<String, Vec<u8>> = load_bam_sequences(bam_path, &opts);
    println!("  {} reads loaded", read_sequences.len());

    println!("Phasing strains");
    let phaser = phase::Phaser::new(bam_path, &target_names, &target_intervals, &read_sequences, output_dir, &opts);
    let haplotypes = phaser.phase(&variants);
    println!("  {} haplotypes phased", haplotypes.len());

    if opts.phase_only {
        println!("Finished!");
        println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
        std::process::exit(0);
    }

    println!("Building aware contigs");
    let aware_contigs = strainberry::awarecontig::build_aware_contigs(&target_sequences, &target_intervals, &haplotypes, opts.min_aware_ctg_len);
    println!("  {} aware contigs built", aware_contigs.len());

    

    println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
}
