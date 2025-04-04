#![allow(dead_code)]

use std::fs;
use std::path::Path;
use std::time::Instant;

use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use anyhow::{bail,Context};
use clap::Parser;
use itertools::Itertools;

use strainberry::bam;
use strainberry::cli;
use strainberry::graph::awaregraph::AwareGraph;
use strainberry::misassembly;
use strainberry::phase;
use strainberry::polish::racon_polish;
use strainberry::seq::build_succinct_sequences;
use strainberry::seq::alignment::load_bam_alignments;
use strainberry::seq::read::{build_read_index,load_bam_sequences};
use strainberry::utils;
use strainberry::variant;


fn main() -> anyhow::Result<(), anyhow::Error> {
    
    let t_start = Instant::now();

    let mut opts = cli::Options::parse();

    utils::check_dependencies(&["minimap2", "samtools", "longcallD", "racon"])?;

    let reads_path = match (opts.in_hifi.as_ref(), opts.in_ont.as_ref()) {
        (None, None) => { bail!("Either --in-hifi or --in-ont is required. For more information, try '--help'.") },
        (Some(hifi_path), None) => {
            opts.mode = cli::Mode::Hifi;
            Path::new(hifi_path)
        },
        (None, Some(ont_path)) => {
            opts.mode = cli::Mode::Nano;
            Path::new(ont_path)
        },
        _ => unreachable!()
    };

    let reference_path = Path::new(&opts.reference);
    // let bam_path = Path::new(&opts.bam);
    let output_dir = Path::new(&opts.output_dir);

    // TODO: consider creating output directory in cli module (after option validation)
    // TODO: consider adding "--force" option and terminate if directory already exists
    fs::create_dir_all(output_dir).with_context(|| format!("Cannot create output directory: \"{}\"", output_dir.display()))?;

    // -------------
    // PREPROCESSING
    // -------------
    
    let preprocess_dir = output_dir.join("00-preprocess");

    println!("Purging duplication from reference: {}", reference_path.display());
    let reference_path = strainberry::derep::derep_assembly(reference_path, &preprocess_dir, &opts)?;
    println!("  dereplicated assembly written to: {}", reference_path.display());

    println!("Mapping reads against the de-replicated assembly");
    let bam_path = preprocess_dir.join("alignment.bam");
    utils::run_minimap2(&reference_path, reads_path, &bam_path, &opts)?;
    println!("  alignment file written to: {}", bam_path.display());

    println!("Loading sequences from: {}", reference_path.display());
    let (target_names, target_index) = bam::bam_target_names(&bam_path);
    let target_sequences = bam::load_sequences(&reference_path, &bam_path);
    println!("  {} sequences loaded", target_sequences.len());

    // let variants = if let Some(vcf) = &opts.vcf {
    //     println!("Loading and filtering variants from: {vcf}");
    //     variant::load_variants_from_vcf(Path::new(vcf), &bam_path, &opts)
    // } else {
    //     println!("Looking for potential variants");
    //     variant::load_variants_from_bam(&bam_path, &opts)
    // };
    let vcf_path = preprocess_dir.join("variants.vcf.gz");
    let variants = variant::run_longcalld(&reference_path, &bam_path, &vcf_path, &opts)?;
    println!("  {} variants found", variants.values().map(|vars| vars.len()).sum::<usize>());

    // TODO:
    // Consider estimating lookback length when a flag "--auto-lookback" is provided
    // utils::estimate_lookback(bam_path, 1000)
    println!("Lookback {} bp", opts.lookback);

    println!("Loading reads from BAM");
    let read_index: HashMap<String, usize> = build_read_index(&bam_path);
    let read_sequences = load_bam_sequences(&bam_path, &read_index, &opts);
    println!("  {} reads loaded", read_index.len());

    let target_intervals = if opts.no_split { 
        bam::bam_target_intervals(&bam_path).into_iter().collect_vec()
    } else {
        println!("Splitting reference at putative misjoins");
        let target_intervals = misassembly::partition_reference(&bam_path, &target_index, &read_index, &opts);
        println!("  {} sequences after split", target_intervals.len());
        target_intervals
    };

    let phased_dir = output_dir.join("20-phased");
    println!("Phasing strains");
    let phaser = phase::Phaser::new(&bam_path, &target_names, &target_intervals, &read_index, &read_sequences, phased_dir, &opts).unwrap();
    let haplotypes = phaser.phase(&variants);
    println!("  {} haplotypes phased", haplotypes.len());

    if opts.phase_only {
        println!("Finished!");
        println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());
        return Ok(());
    }

    println!("Building strain-aware contigs");
    let mut aware_contigs = strainberry::awarecontig::build_aware_contigs(&target_intervals, &haplotypes, opts.min_aware_ctg_len);
    println!("  {} strain-aware contigs built", aware_contigs.len());

    println!("Loading read alignments");
    let read_alignments = load_bam_alignments(&bam_path, &read_index, &opts);

    println!("Building succinct reads");
    let succinct_reads = build_succinct_sequences(&bam_path, &variants, &read_index, &opts);

    println!("Read realignment to haplotypes");
    let variant_positions = variants.values().flatten().map(|var| (var.tid,var.pos)).collect::<HashSet<_>>();
    let (seq2haplo, ambiguous_reads) = strainberry::phase::separate_reads(&succinct_reads, &haplotypes, &variant_positions, opts.min_shared_snv);
    println!("  unique={} ambiguous={}", seq2haplo.len(), ambiguous_reads.len());

    println!("Mapping reads to strain-aware contigs");
    let read2aware = strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut aware_contigs, &seq2haplo, &ambiguous_reads);

    let graphs_dir = output_dir.join("40-graphs");
    fs::create_dir_all(graphs_dir.as_path()).with_context(|| format!("Cannot create graphs directory: \"{}\"", graphs_dir.display()))?;
    
    println!("Building strain-aware graph");
    let mut aware_graph = AwareGraph::build(&aware_contigs);
    aware_graph.add_edges_from_aware_alignments(&read2aware);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.raw.gfa"), &target_names).unwrap();
    aware_graph.write_dot(graphs_dir.join("aware_graph.raw.dot")).unwrap();

    aware_graph.remove_weak_edges(5);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.gfa"), &target_names).unwrap();

    println!("Resolving strain-aware graph using reads");
    let nb_tedges = aware_graph.add_bridges(&read2aware);
    aware_graph.write_dot(graphs_dir.join("aware_graph.dot")).unwrap();
    println!("  {} read bridges added", nb_tedges);

    let mut num_iter = 1;
    let mut num_resolved = aware_graph.resolve_read_bridges(opts.min_alt_count);
    let mut tot_resolved = num_resolved;
    while num_resolved > 0 {
        num_iter += 1;
        num_resolved = aware_graph.resolve_read_bridges(opts.min_alt_count);
        tot_resolved += num_resolved;
    }
    println!("  {tot_resolved} junctions resolved after {num_iter} iterations");
    aware_graph.write_gfa(graphs_dir.join("aware_graph.resolved.gfa"), &target_names)?;
    aware_graph.clear_transitive_edges();

    println!("Building assembly graph");
    let unitig_dir = output_dir.join("50-unitigs");
    let unitig_graph = aware_graph.build_assembly_graph(&target_names, &target_sequences, &read_sequences, phaser.fragments_dir(), &unitig_dir, &opts)?;
    unitig_graph.write_gfa(&graphs_dir.join("final_unpolished.gfa"))?;

    println!("Polishing");
    let polish_dir = output_dir.join("60-polishing");
    let unpolished_fasta_path = graphs_dir.join("final_unpolished.fasta");
    unitig_graph.write_fasta(&unpolished_fasta_path)?;
    let in_reads_path = Path::new(opts.in_hifi.as_ref().unwrap());
    let polished_fasta_path = output_dir.join("assembly.fasta");
    racon_polish(&unpolished_fasta_path, in_reads_path, &polished_fasta_path, strainberry::polish::PolishMode::Aware, &polish_dir, &opts)?;

    println!("  polished assembly written to {}", polished_fasta_path.display());

    println!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());

    Ok(())
}
