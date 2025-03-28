use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;


pub fn purge_dups(contigs: &Path, output: &Path) -> Result<()> {
    // let file = match File::open(path) {
    //     Err(err) => panic!("Could not open file \"{}\": {}", path.display(), err),
    //     Ok(file) => file,
    // };
    // match path.extension() {
    //     Some(ext) if ext == "gz" => Box::new(BufReader::new(MultiGzDecoder::new(file))),
    //     _ => Box::new(BufReader::new(file)),
    // }
    Ok(())
}