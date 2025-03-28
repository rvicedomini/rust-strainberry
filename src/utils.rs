use std::fs::File;
use std::io::{BufRead,BufReader,Write,BufWriter};
use std::mem::MaybeUninit;
use std::path::Path;

use ahash::AHashMap as HashMap;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

pub fn get_maxrss() -> f64 {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as f64 / (1024.0 * 1024.0)
}


pub fn get_file_reader(path: &Path) -> Box<dyn BufRead> {
    let file = match File::open(path) {
        Err(err) => panic!("Could not open file \"{}\": {}", path.display(), err),
        Ok(file) => file,
    };
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    }
}


pub fn get_file_writer(path: &Path) -> Box<dyn Write> {
    let file = match File::create(path) {
        Err(err) => panic!("Could not open/create file \"{}\": {}", path.display(), err),
        Ok(file) => file,
    };
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(
            BufWriter::new(GzEncoder::new(file, Compression::default()))
        ),
        _ => Box::new(BufWriter::new(file)),
    }
}


// flip strand as u8 character
// b'+': 00101|01|1
// b'-': 00101|10|1
//    6: 00000|11|0
#[inline(always)]
pub fn flip_strand(strand:u8) -> u8 {
    assert!(strand == b'+' || strand == b'-');
    strand ^ 6
}

#[inline(always)]
pub fn counter_from_iter<T>(iter: impl Iterator<Item = T>) -> HashMap<T,usize> 
where
    T: std::cmp::Eq + std::hash::Hash
{    
    let mut counter = HashMap::new();
    for e in iter {
        counter.entry(e).and_modify(|e| *e += 1).or_insert(1);
    }
    counter
}
