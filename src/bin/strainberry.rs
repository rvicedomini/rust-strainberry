#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;

// use ahash::AHashMap as HashMap;
use anyhow::{bail,Context};
use clap::Parser;
use itertools::Itertools;

use strainberry::cli;
use strainberry::graph::awaregraph::AwareGraph;
use strainberry::misassembly;
use strainberry::phase;
// use strainberry::polish::racon_polish;
use strainberry::seq::{self, SeqInterval};
use strainberry::utils;
use strainberry::variant;


fn main() -> anyhow::Result<(), anyhow::Error> {

    spdlog::default_logger().set_level_filter(spdlog::LevelFilter::All);
    
    let t_start = Instant::now();

    let mut opts = cli::Options::parse();

    utils::check_dependencies(&["minimap2", "samtools", "racon"])?;

    if opts.caller == cli::VarCaller::Longcalld {
        utils::check_dependencies(&["longcallD"])?;
    }

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
        spdlog::info!("Purging duplications from reference: {}", reference_path.display());
        reference_path = strainberry::derep::derep_assembly(&reference_path, &preprocess_dir, &opts)?;
    }

    let bam_path = if !opts.no_derep || opts.bam.is_none() {
        spdlog::info!("Mapping reads to reference assembly: {}", reference_path.display());
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
    let read_db = seq::SeqDatabase::build(reads_path, false)?;
    spdlog::info!("{} sequences processed", read_db.size());

    let variants = if let Some(vcf) = &opts.vcf {
        spdlog::info!("Loading and filtering variants from: {vcf}");
        variant::load_variants_from_vcf(Path::new(vcf), &bam_path, &ref_db, &opts)
    } else if opts.caller == cli::VarCaller::Longcalld {
        spdlog::info!("Calling variants using longcallD");
        let vcf_path = preprocess_dir.join("variants.vcf.gz");
        variant::run_longcalld(&reference_path, &bam_path, &ref_db, &vcf_path, &opts)?
    } else {
        spdlog::info!("Calling variants from pileup");
        variant::load_variants_from_bam(&bam_path, &ref_db, &opts)
    };
    spdlog::info!("{} variants found", variants.values().map(|vars| vars.len()).sum::<usize>());    

    // TODO:
    // Consider setting lookback as read Nx for an appropriate value of x.
    let read_n75 = strainberry::bam::estimate_lookback(&bam_path, 75, 0).unwrap();
    spdlog::info!("Read N75 {} bp", read_n75);
    // opts.lookback = read_n75;

    let ref_intervals = if opts.no_split {
        ref_db.sequences.iter().enumerate()
            .map(|(ref_idx,seq)| SeqInterval{ tid:ref_idx, beg:0, end:seq.len() })
            .collect_vec()
    } else {
        spdlog::info!("Splitting reference at putative misjoins");
        let intervals = misassembly::partition_reference(&bam_path, &ref_db, &read_db, &opts);
        spdlog::info!("{} sequences after split", intervals.len());
        intervals
    };

    spdlog::info!("Loading read alignments");
    let read_alignments = seq::alignment::load_bam_alignments(&bam_path, &ref_db, &read_db, &opts);

    // if opts.no_phase {

    //     let mut aware_contigs = strainberry::awarecontig::build_aware_contigs(&ref_intervals, &HashMap::new(), opts.min_aware_ctg_len);
        
    //     spdlog::info!("{} strain-aware contigs built", aware_contigs.len());

    //     spdlog::info!("Building succinct reads");
    //     let succinct_reads = seq::build_succinct_sequences(&bam_path, &ref_db, &read_db, &variants, &opts);

    //     spdlog::info!("Read realignment to haplotypes");
    //     let seq2haplo = strainberry::phase::separate_reads(&succinct_reads, &HashMap::new(), opts.min_shared_snv);

    //     spdlog::info!("Mapping reads to strain-aware contigs");
    //     let read2aware = strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut aware_contigs, &seq2haplo);

    //     let graphs_dir = output_dir.join("40-graphs");
    //     fs::create_dir_all(graphs_dir.as_path()).with_context(|| format!("Cannot create graphs directory: \"{}\"", graphs_dir.display()))?;
        
    //     spdlog::info!("Building strain-aware graph");
    //     let mut aware_graph = AwareGraph::build(&aware_contigs);
    //     aware_graph.add_edges_from_aware_alignments(&read2aware);
    //     aware_graph.write_gfa(graphs_dir.join("aware_graph.raw.gfa"), &ref_db)?;

    //     aware_graph.remove_weak_edges(5);
    //     aware_graph.write_gfa(graphs_dir.join("aware_graph.gfa"), &ref_db)?;

    //     spdlog::info!("Strain-aware graph resolution");
    //     let nb_tedges = aware_graph.add_bridges(&read2aware);
    //     aware_graph.write_dot(graphs_dir.join("aware_graph.dot"))?;
    //     spdlog::debug!("{nb_tedges} read bridges added");

    //     let nb_resolved = aware_graph.resolve_read_bridges(opts.min_alt_count);
    //     spdlog::info!("{nb_resolved} junctions resolved");

    //     aware_graph.write_gfa(graphs_dir.join("aware_graph.resolved.gfa"), &ref_db)?;
    //     aware_graph.write_dot(graphs_dir.join("aware_graph.resolved.dot"))?;

    //     spdlog::info!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
    //     return Ok(());
    // }

    let phased_dir = output_dir.join("20-phased");
    spdlog::info!("Phasing strains");
    let phaser = phase::Phaser::new(&bam_path, &ref_db, &read_db, &ref_intervals, phased_dir, &opts).unwrap();
    let haplotypes = phaser.phase(&variants);
    spdlog::info!("{} haplotypes phased", haplotypes.len());

    spdlog::info!("Building strain-aware contigs");
    let mut aware_contigs = strainberry::awarecontig::build_aware_contigs(&ref_intervals, &haplotypes, opts.min_aware_ctg_len);
    spdlog::info!("{} strain-aware contigs built", aware_contigs.len());

    spdlog::info!("Building succinct reads");
    let succinct_reads = seq::build_succinct_sequences(&bam_path, &ref_db, &read_db, &variants, &opts);

    spdlog::info!("Read realignment to haplotypes");
    let seq2haplo = strainberry::phase::separate_reads(&succinct_reads, &haplotypes, opts.min_shared_snv);

    spdlog::info!("Mapping reads to strain-aware contigs");
    let read2aware = strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut aware_contigs, &seq2haplo);

    let graphs_dir = output_dir.join("40-graphs");
    fs::create_dir_all(graphs_dir.as_path()).with_context(|| format!("Cannot create graphs directory: \"{}\"", graphs_dir.display()))?;
    
    spdlog::info!("Building strain-aware graph");
    let mut aware_graph = AwareGraph::build(&aware_contigs);
    aware_graph.add_edges_from_aware_alignments(&read2aware);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.raw.gfa"), &ref_db)?;
    aware_graph.write_dot(graphs_dir.join("aware_graph.raw.dot"))?;

    aware_graph.remove_weak_edges(opts.min_alt_count);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.gfa"), &ref_db)?;

    spdlog::info!("Strain-aware graph resolution");
    let nb_tedges = aware_graph.add_bridges(&read2aware);
    aware_graph.write_dot(graphs_dir.join("aware_graph.dot"))?;
    spdlog::debug!("{nb_tedges} read bridges added");

    let nb_resolved = aware_graph.resolve_read_bridges(opts.min_alt_count);
    spdlog::info!("{nb_resolved} junctions resolved");

    aware_graph.write_gfa(graphs_dir.join("aware_graph.resolved.gfa"), &ref_db)?;
    aware_graph.write_dot(graphs_dir.join("aware_graph.resolved.dot"))?;
    aware_graph.clear_transitive_edges();

    if opts.no_asm {
        spdlog::info!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
        return Ok(());
    }

    spdlog::info!("Building assembly graph");
    let unitig_dir = output_dir.join("50-unitigs");
    let unitig_graph = aware_graph.build_assembly_graph(&ref_db, &read_db, phaser.fragments_dir(), &unitig_dir, &opts)?;
    unitig_graph.write_gfa(&output_dir.join("assembly.gfa"))?;
    let assembly_fasta_path = output_dir.join("assembly.fasta");
    unitig_graph.write_fasta(&assembly_fasta_path)?;
    spdlog::info!("Final assembly written to {}", assembly_fasta_path.display());

    // spdlog::info!("Polishing assembly");
    // let polish_dir = output_dir.join("60-polishing");
    // let in_reads_path = Path::new(opts.in_hifi.as_ref().unwrap());
    // let polished_fasta_path = output_dir.join("assembly.fasta");
    // racon_polish(&unpolished_fasta_path, in_reads_path, &polished_fasta_path, strainberry::polish::PolishMode::Aware, &polish_dir, &opts)?;
    // spdlog::info!("Final assembly written to {}", polished_fasta_path.display());

    spdlog::info!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());

    Ok(())
}
