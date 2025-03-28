// use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::str::FromStr;

use anyhow::{Context, Error, Result};
use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use itertools::Itertools;
// use flate2::read::MultiGzDecoder;
// use itertools::Itertools;
use needletail::Sequence;

use crate::cli::Options;
use crate::seq::alignment::{MappingType, PafAlignment};


fn map_type(a: &PafAlignment, overhang: usize, r: f64) -> MappingType {
    let query_range = if a.strand == b'+' {
        (a.query_beg, a.query_end, a.query_length)
    } else {
        (a.query_length-a.query_end, a.query_length-a.query_beg, a.query_length)
    };
    let target_range = (a.target_beg, a.target_end, a.target_length);
    crate::seq::alignment::classify_mapping(query_range, target_range, overhang, r)
}

fn compute_matching_intervals(fasta_path: &Path, opts: &Options) -> Result<HashMap<String,Vec<(usize,usize,String)>>,Error> {

    let fasta_path_str = fasta_path.to_str().unwrap();

    // if which::which("minimap2").is_err() {
    //     bail!("cannot execute minimap2, please check it is in your system PATH");
    // }

    let args = ["-t", &opts.nb_threads.to_string(), "-cDP", "--dual=no", "--no-long-join", "-r85", fasta_path_str, fasta_path_str];
    let stdout = Command::new("minimap2").args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn().context("cannot run minimap2")?
        .stdout
        .with_context(|| "Could not capture standard output.")?;
    
    let reader = BufReader::new(stdout);
    let alignments: Vec<PafAlignment> = reader.lines()
        .map_while(Result::ok)
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty())
        .filter_map(|line| PafAlignment::from_str(&line).ok())
        .collect();

    let mut matching_intervals: HashMap<String,Vec<(usize,usize,String)>> = HashMap::new();
    for a in alignments {
        let map_type = self::map_type(&a, 200, 0.5);
        let identity = if a.mapping_length > 0 { (100.0 * a.matches as f64)/(a.mapping_length as f64) } else { 0.0 };
        if a.query_name == a.target_name || identity < 95.0 {
            continue
        }
        match map_type {
            MappingType::QueryContained => {
                matching_intervals.entry(a.query_name).or_default().push((0,a.query_length,a.target_name));
            },
            MappingType::QueryPrefix => {
                let (query_beg, query_end) = if a.strand == b'+' { (0,a.query_end) } else { (a.query_beg,a.query_length) };
                matching_intervals.entry(a.query_name).or_default().push((query_beg, query_end,a.target_name));
            },
            MappingType::QuerySuffix => {
                let (query_beg, query_end) = if a.strand == b'+' { (a.query_beg,a.query_length) } else { (0,a.query_end) };
                matching_intervals.entry(a.query_name).or_default().push((query_beg, query_end,a.target_name));
            },
            MappingType::ReferenceContained => {
                matching_intervals.entry(a.target_name).or_default().push((0,a.target_length,a.query_name));
            },
            MappingType::ReferencePrefix => {
                matching_intervals.entry(a.target_name).or_default().push((0,a.target_end,a.query_name));
            },
            MappingType::ReferenceSuffix => {
                matching_intervals.entry(a.target_name).or_default().push((a.target_beg,a.target_length,a.query_name));
            },
            MappingType::DovetailPrefix => {
                let (query_beg, query_end) = if a.strand == b'+' { (0,a.query_end) } else { (a.query_beg,a.query_length) };
                matching_intervals.entry(a.query_name.clone()).or_default().push((query_beg, query_end,a.target_name.clone()));
                matching_intervals.entry(a.target_name).or_default().push((a.target_beg, a.target_length, a.query_name));
            },
            MappingType::DovetailSuffix => {
                let (query_beg, query_end) = if a.strand == b'+' { (a.query_beg,a.query_length) } else { (0,a.query_end) };
                matching_intervals.entry(a.query_name.clone()).or_default().push((query_beg, query_end,a.target_name.clone()));
                matching_intervals.entry(a.target_name).or_default().push((0, a.target_end, a.query_name));
            },
            MappingType::Internal => {}
        }
    }

    // let mut fasta_reader = needletail::parse_fastx_file(fasta_path)
    //     .with_context(|| format!("Cannot open fasta file: {}", fasta_path.display()))?;

    Ok(matching_intervals)
}


pub fn derep_assembly(fasta_path: &Path, work_dir: PathBuf, opts: &Options) -> Result<PathBuf> {
    assert!(fasta_path.is_file());

    let tmp_dir = work_dir.join("tmp");
    std::fs::create_dir_all(&tmp_dir)
        .with_context(|| format!("Cannot create output directory: \"{}\"", tmp_dir.display()))?;

    let output_path = work_dir.join("assembly_derep.fasta");
    let mut asm_path = fasta_path.to_path_buf();
    
    for it in 1.. {

        let mut modified = false;
        let matching_intervals = self::compute_matching_intervals(&asm_path, opts)?;

        eprintln!("minimap2 successfully run in iteration {it}");

        let mut seq_dict: HashMap<String,String> = HashMap::new();
        let mut fasta_reader = needletail::parse_fastx_file(fasta_path).context("Cannot open fasta file")?;
        while let Some(record) = fasta_reader.next() {
            let record = record.unwrap();
            let name = std::str::from_utf8(record.id()).unwrap().to_string();
            let sequence = std::str::from_utf8(&record.normalize(false)).unwrap().to_string();
            seq_dict.insert(name, sequence);
        }

        let derep_path = tmp_dir.join(format!("derep_{it}.fasta"));
        let mut derep_writer = crate::utils::get_file_writer(&derep_path);
        let mut processed: HashSet<&str> = HashSet::new();
        for (seq_id, sequence) in seq_dict.iter().sorted_unstable_by_key(|(_,seq)| seq.len()) {
            let mut seq_intervals = vec![(0,sequence.len())];
            if let Some(intervals) = matching_intervals.get(seq_id.as_str()) {
                let sorted_intervals = intervals.iter()
                    .filter_map(|(beg,end,name)| {
                        if processed.contains(name.as_str()) { Some((*beg,*end,name.clone())) } else { None }
                    }).sorted_unstable().collect_vec();
                
                let mut merged_intervals = Vec::with_capacity(sorted_intervals.len());
                for iv in sorted_intervals {
                    if merged_intervals.is_empty() {
                        merged_intervals.push(iv);
                        continue;
                    }
                    let last = unsafe { merged_intervals.last_mut().unwrap_unchecked() };
                    if iv.0 <= last.1 { last.1 = std::cmp::max(last.1, iv.1); }
                }
                modified = !merged_intervals.is_empty();
                for iv in &merged_intervals {
                    assert!(iv.0 < iv.1);
                    let first = seq_intervals.partition_point(|(_,end)| *end <= iv.0);
                    let last = seq_intervals.partition_point(|(beg,_)| *beg < iv.1 );
                    let overlapping = seq_intervals.drain(first..last).collect_vec();
                    for (beg,end) in overlapping {
                        if iv.0 <= beg && end <= iv.1 {
                            continue
                        } else if beg <= iv.0 && iv.1 <= end {
                            if beg < iv.0 { seq_intervals.push((beg,iv.0)); }
                            if iv.1 < end { seq_intervals.push((iv.1,end)); }
                        } else if beg < iv.0 {
                            assert!(end < iv.1);
                            seq_intervals.push((beg,iv.0))
                        } else {
                            assert!(iv.0 < beg && iv.1 < end);
                            seq_intervals.push((iv.1,end));
                        }
                    }
                }
                processed.insert(seq_id.as_str());
            }

            // Sequences are processed in increasing length and only alignments involving the suffix/prefix are retained
            // This means I cannot get more than 1 interval
            assert!(seq_intervals.len() <= 1);
            for (beg,end) in seq_intervals {
                let record = format!(">{seq_id}\n{}\n", &sequence[beg..end]);
                derep_writer.write_all(record.as_bytes())?;
            }
        }

        if !modified {
            std::fs::rename(&derep_path, &output_path)
                .with_context(|| format!("cannot move file \"{}\" to \"{}\"", derep_path.display(), output_path.display()))?;
            break
        }

        asm_path = derep_path;
    }

    Ok(output_path)
}

// def main(argv=None):
    
//     with TemporaryDirectory(dir='.', prefix='derep') as tmpdir:

//         logger.debug(f'tmpdir: {tmpdir}')
    
//         asmfile = opt.ASM
//         modified = True
//         while modified:

//             modified = False
//             with derep_asmfile as derep:
//                 processed = set()
//                 for seq_id, sequence in sorted(seq_dict.items(), key=lambda x:len(x[1])):
//                     seq_length = len(sequence)
//                     seq_intervals = IntervalTree([Interval(0,seq_length)])
//                     if seq_id in matching_intervals:
//                         matching_intervals[seq_id] = IntervalTree.from_tuples((iv.begin,iv.end,iv.data) for iv in matching_intervals[seq_id] if (iv.data not in processed))
//                         matching_intervals[seq_id].merge_overlaps()
//                         for iv in matching_intervals[seq_id]:
//                             assert(iv.end > iv.begin)
//                             seq_intervals.chop(iv.begin,iv.end)
//                             modified = True
//                         processed.add(seq_id)
//                     assert(len(seq_intervals) <= 1)
//                     for iv in seq_intervals:
//                         if iv.end - iv.begin >= opt.min_length:
//                             derep.write(f'>{seq_id}\n{insert_newlines(sequence[iv.begin:iv.end])}\n')
            
//             if not modified:
//                 Path(derep_asmfile.name).rename(opt.outfile)
//                 logger.debug(f'outfile: {opt.outfile}')
//                 break

//             asmfile = derep_asmfile.name
    
//     return 0