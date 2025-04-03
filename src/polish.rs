use std::io::{BufRead,BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::str::FromStr;

use anyhow::{bail, Context, Result};
use ahash::AHashMap as HashMap;

use crate::cli::Options;
use crate::seq::alignment::{MappingType, PafAlignment};


fn filter_paf_alignments() {}

pub fn racon_polish(target_path: &Path, read_path: &Path, work_dir: PathBuf, opts: &Options) -> Result<PathBuf> {

    if which::which("racon").is_err() {
        bail!("missing racon dependency, please check your system PATH");
    }

    std::fs::create_dir_all(&work_dir)
        .with_context(|| format!("cannot create output directory: \"{}\"", work_dir.display()))?;

    // Run minimap2 and load alignments

    let mm2_preset = match opts.mode {
        crate::cli::Mode::Hifi => "-xmap-hifi",
        crate::cli::Mode::Nano => "-xmap-ont",
    };
    
    let args = ["-t", &opts.nb_threads.to_string(), "-c", mm2_preset, target_path.to_str().unwrap(), read_path.to_str().unwrap()];
    let stdout = Command::new("minimap2").args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn().context("cannot run minimap2")?
        .stdout
        .with_context(|| "Could not capture standard output.")?;
    
    let mut alignments: HashMap<String,PafAlignment> = HashMap::new();
    let reader = BufReader::new(stdout);
    reader.lines()
        .map_while(Result::ok)
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty())
        .map(|line| PafAlignment::from_str(&line).with_context(||format!("error parsing line:\n{line}")).unwrap())
        .filter(|paf| matches!(paf.map_type(10, 0.01), MappingType::DovetailPrefix|MappingType::DovetailSuffix|MappingType::QueryContained|MappingType::ReferenceContained))
        .for_each(|paf| {
            let qname = paf.query_name.clone();
            if !alignments.contains_key(&qname) {
                alignments.insert(qname, paf);
            } else {
                let e = alignments.get_mut(&qname).unwrap();
                if paf.identity > e.identity { *e = paf; }
            }
        });

    // Write retained alignments in PAF format

    let paf_path = {
        let paf_path = work_dir.join("alignments.paf.gz");
        let mut paf_writer = crate::utils::get_file_writer(&paf_path);
        for paf in alignments.into_values() {
            paf_writer.write_all(paf.to_string().as_bytes())?;
            paf_writer.write_all(b"\n")?;
        }
        paf_path
    };

    // Run racon

    let polished_path = work_dir.join("polished.fasta");
    let polished_file = std::fs::File::create(&polished_path).unwrap();
    let racon_log_path = work_dir.join("racon.log");
    let racon_log_file = std::fs::File::create(&racon_log_path).unwrap();

    let args = ["-t", &opts.nb_threads.to_string(), "-u", read_path.to_str().unwrap(), paf_path.to_str().unwrap(), target_path.to_str().unwrap()];
    // println!("running cmd: racon {}", args.join(" "));
    let mut racon = Command::new("racon").args(args)
        .stdout(Stdio::from(polished_file))
        .stderr(Stdio::from(racon_log_file))
        .spawn().context("cannot run racon")?;
    
    racon.wait().unwrap();

    Ok(polished_path)
}
