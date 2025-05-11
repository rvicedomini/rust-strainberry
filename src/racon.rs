mod alignment;
mod polish;
mod sequence;
mod window;

// use std::path::Path;
// use std::time::Instant;


// fn main() -> anyhow::Result<(), anyhow::Error> {

//     spdlog::default_logger().set_level_filter(spdlog::LevelFilter::All);
    
//     let t_start = Instant::now();

//     let opts = cli::Options::parse();

//     let reference_path = Path::new(&opts.reference);
//     let reads_path = Path::new(&opts.reads);
//     let paf_path = Path::new(&opts.paf);

//     rayon::ThreadPoolBuilder::new().num_threads(opts.nb_threads).build_global().unwrap();

//     let (ref_sequences, ref_index) = sequence::load_sequences(reference_path);
//     info!("loaded reference sequences: {}", ref_sequences.len());

//     let (read_sequences, read_index) = sequence::load_sequences(reads_path);
//     info!("loaded read sequences: {}", read_sequences.len());

//     let alignments = alignment::load_paf_alignments(paf_path, &ref_index, &read_index)?;
//     let alignments = alignment::filter_alignments(alignments, &opts);
//     if alignments.is_empty() {
//         warn!("empty overlap set, no output file generated");
//         return Ok(())
//     }
//     info!("loaded alignments: {}", alignments.len());

//     let mut polisher = polish::Polisher::new(&ref_sequences, &read_sequences, alignments, &opts);

//     info!("initializing polisher");
//     polisher.initialize();

//     info!("polishing windows");
//     let polished_sequences = polisher.polish();

//     let out_path = Path::new(&opts.output);
//     let mut writer = crate::utils::get_file_writer(out_path);
//     for (i, seq) in polished_sequences.into_iter().enumerate() {
//         let header = format!(">{i}\n");
//         writer.write_all(header.as_bytes())?;
//         let sequence = utils::insert_newlines(std::str::from_utf8(&seq).unwrap(),120);
//         writer.write_all(sequence.as_bytes())?;
//         writer.write_all(b"\n")?;
//     }
    
//     info!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());

//     Ok(())
// }
