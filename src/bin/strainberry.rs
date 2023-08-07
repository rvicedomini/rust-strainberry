use std::fs;
use std::path::Path;
use std::time::Instant;

use clap::Parser;

use strainberry::cli;
use strainberry::misassembly;
use strainberry::phase;
use strainberry::utils;
use strainberry::variant;


fn main() {
    
    let t_start = Instant::now();

    let opts = cli::Options::parse();

    let fasta_path = Path::new(&opts.fasta_file);
    let bam_path = Path::new(&opts.bam_file);
    let output_dir = Path::new(&opts.output_dir);

    // TODO?
    // create output directory in cli module, after option validation
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
    println!("Loaded {} variants", variants.values().map(|vars| vars.len()).sum::<usize>());

    println!("Loading sequences from: {}", fasta_path.canonicalize().unwrap().display());
    let target_sequences: Vec<Vec<u8>> = utils::load_sequences(fasta_path, bam_path);
    println!("Loaded {} sequences", target_sequences.len());
    
    // TODO:
    // Consider estimating lookback length when a flag "--auto-lookback" is provided
    // utils::estimate_lookback(bam_path, 1000)
    println!("Lookback {} bp", opts.lookback);

    misassembly::find_misassemblies(bam_path, &opts);

    // println!("Phasing strains");
    // let phaser = phase::Phaser::new(&bam_path, &target_sequences, &output_dir, &opts);
    // phaser.run(&variants);

    // separate_workdir = os.path.join(opt.outdir,'20-separate')
    // separator = HifiReadSeparator(reference_lengths, opt.BAM, read_dict, variant_positions, reference_alignments, reference_intervals, separate_workdir,
    //     mapq=opt.min_mapq, lookback=opt.lookback, min_obs=opt.min_read_obs, min_frac=opt.min_read_frac,
    //     graph_only=opt.graph_only, phase_only=opt.phase_only, debug=opt.debug)
    // separator.separate_references()

    println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
}
