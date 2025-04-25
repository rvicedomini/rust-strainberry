#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;

use ahash::AHashMap as HashMap;
use anyhow::{bail,Context};
use clap::Parser;
use itertools::Itertools;

use strainberry::cli;
use strainberry::graph::awaregraph::AwareGraph;
use strainberry::misassembly;
use strainberry::phase;
use strainberry::seq;
use strainberry::utils;
use strainberry::variant;


fn main() -> anyhow::Result<(), anyhow::Error> {
    
    let t_start = Instant::now();

    let opts = cli::Options::parse();

    if opts.debug {
        spdlog::default_logger().set_level_filter(spdlog::LevelFilter::MoreSevereEqual(spdlog::Level::Debug));
    }

    utils::check_dependencies(&["minimap2", "samtools", "racon"])?;
    // if opts.caller == cli::VarCaller::Longcalld {
    //     utils::check_dependencies(&["longcallD"])?;
    // }

    if opts.in_hifi.is_none() && opts.in_ont.is_none() {
        bail!("Either --in-hifi or --in-ont is required. For more information, try '--help'.")
    }

    run_pipeline(opts)?;

    spdlog::info!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());

    Ok(())
}


fn run_pipeline(mut opts: cli::Options) -> anyhow::Result<(), anyhow::Error> {
    
    let reads_path = match (opts.in_hifi.as_ref(), opts.in_ont.as_ref()) {
        (None, None) => { bail!("Either --in-hifi or --in-ont is required. For more information, try '--help'.") },
        (Some(hifi_path), None) => {
            opts.mode = cli::Mode::Hifi;
            opts.strand_bias_pvalue = 0.005;
            Path::new(hifi_path)
        },
        (None, Some(ont_path)) => {
            opts.mode = cli::Mode::Nano;
            opts.strand_bias_pvalue = 0.01;
            Path::new(ont_path)
        },
        _ => unreachable!()
    };

    let mut reference_path = PathBuf::from_str(&opts.reference)?;
    let output_dir = Path::new(&opts.output_dir);

    // TODO: consider creating output directory in cli module (after option validation)
    // TODO: consider adding "--force" option and terminate if directory already exists
    fs::create_dir_all(output_dir).with_context(|| format!("Cannot create output directory: \"{}\"", output_dir.display()))?;

    // -------------
    // PREPROCESSING
    // -------------
    
    let preprocess_dir = output_dir.join("00-preprocess");
    std::fs::create_dir_all(&preprocess_dir)
        .with_context(|| format!("Cannot create output directory: \"{}\"", preprocess_dir.display()))?;

    if let Some("gfa") = reference_path.extension().and_then(std::ffi::OsStr::to_str) {
        spdlog::info!("Converting input GFA to FASTA");
        let gfa_graph = strainberry::graph::GfaGraph::from_gfa(&reference_path)?;
        let mut fasta_name = reference_path.file_stem().unwrap().to_os_string();
        fasta_name.push(".fasta");
        reference_path = preprocess_dir.join(fasta_name);
        gfa_graph.write_fasta(&reference_path)?;
    }

    if !opts.no_derep {
        spdlog::info!("Purging duplications from input reference");
        reference_path = strainberry::derep::derep_assembly(&reference_path, &preprocess_dir, &opts)?;
        spdlog::info!("Dereplicated reference written to {}", reference_path.display());
    }

    let bam_path = if !opts.no_derep || opts.bam.is_none() {
        spdlog::info!("Mapping reads to dereplicated reference");
        let bam_path = preprocess_dir.join("alignment.bam");
        utils::run_minimap2(&reference_path, reads_path, &bam_path, &opts)?;
        bam_path
    } else {
        PathBuf::from_str(opts.bam.as_ref().unwrap())?
    };

    spdlog::info!("Loading sequences from: {}", reference_path.display());
    let ref_db = seq::SeqDatabase::build(&reference_path, true)?;    
    spdlog::info!("{} sequences processed", ref_db.size());

    spdlog::info!("Building read index");
    let read_db = seq::SeqDatabase::build(reads_path, true)?;
    let read_n75 = read_db.compute_nx(75).unwrap();
    opts.lookback = (9 * read_n75) / 10;
    spdlog::debug!("{} sequences loaded / N75: {read_n75}", read_db.size());

    spdlog::info!("Calling variants from pileup");
    let variants = variant::load_variants_from_bam(&bam_path, &ref_db, &opts);
    spdlog::debug!("{} variants identified", variants.values().map(|vars| vars.len()).sum::<usize>());

    spdlog::info!("Splitting reference at putative misjoins");
    let ref_intervals = misassembly::partition_reference(&bam_path, &ref_db, &read_db, &opts);
    spdlog::debug!("{} sequences after split", ref_intervals.len());

    spdlog::info!("Loading read alignments");
    let read_alignments = seq::alignment::load_bam_alignments(&bam_path, &ref_db, &read_db, &opts);

    if opts.no_phase {
        return Ok(());
    }

    let ref_intervals = {
        spdlog::info!("Filtering low-coverage sequences");
        let mut ref_contigs = strainberry::awarecontig::build_aware_contigs(&ref_intervals, &HashMap::new(), opts.min_aware_ctg_len);
        strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut ref_contigs, &HashMap::new());
        ref_contigs.retain(|ctg| ctg.depth() >= opts.min_alt_count as f64);
        ref_contigs.into_iter().map(|ctg| ctg.interval()).collect_vec()
    };

    let phased_dir = output_dir.join("20-phased");
    spdlog::info!("Phasing strain haplotypes");
    let phaser = phase::Phaser::new(&bam_path, &ref_db, &read_db, &ref_intervals, phased_dir, &opts).unwrap();
    let phaser_result = phaser.phase(&variants);
    spdlog::info!("{} haplotypes phased", phaser_result.haplotypes.len());

    spdlog::info!("Building strain-aware contigs");
    let mut aware_contigs = strainberry::awarecontig::build_aware_contigs(&ref_intervals, &phaser_result.haplotypes, 0);

    spdlog::info!("Building succinct reads");
    let succinct_reads = seq::build_succinct_sequences(&bam_path, &ref_db, &read_db, &variants, &opts);

    spdlog::info!("Read realignment to haplotypes");
    let seq2haplo = strainberry::phase::separate_reads(&succinct_reads, &phaser_result.haplotypes, opts.min_shared_snv);
    drop(phaser_result);

    spdlog::info!("Mapping reads to strain-aware contigs");
    let read2aware = strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut aware_contigs, &seq2haplo);

    let graphs_dir = output_dir.join("40-graphs");
    fs::create_dir_all(graphs_dir.as_path()).with_context(|| format!("Cannot create graphs directory: \"{}\"", graphs_dir.display()))?;
    
    spdlog::info!("Building strain-aware graph");
    let mut aware_graph = AwareGraph::build(&aware_contigs);
    aware_graph.add_edges_from_aware_alignments(&read2aware);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.raw.gfa"), &ref_db)?;

    spdlog::info!("Strain-aware graph resolution");
    aware_graph.add_bridges(&read2aware);
    aware_graph.remove_weak_edges(opts.min_alt_count);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.gfa"), &ref_db)?;
    
    aware_graph.resolve_junctions(opts.min_alt_count);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.resolved.gfa"), &ref_db)?;

    if opts.no_asm {
        return Ok(());
    }

    spdlog::info!("Building assembly graph");
    let unitig_dir = output_dir.join("50-unitigs");
    let unitig_graph = aware_graph.build_assembly_graph(&ref_db, &read_db, phaser.fragments_dir(), &unitig_dir, &opts)?;
    unitig_graph.write_gfa(&output_dir.join("assembly.gfa"))?;
    let assembly_fasta_path = output_dir.join("assembly.fasta");
    unitig_graph.write_fasta(&assembly_fasta_path)?;
    spdlog::info!("Final assembly written to {}", assembly_fasta_path.display());
    
    Ok(())
}
