use std::mem::MaybeUninit;
use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;
use flate2::read::MultiGzDecoder;


pub fn get_maxrss() -> f64 {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as f64 / (1024.0 * 1024.0)
}


pub fn get_file_reader(path: &Path) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap();
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    }
}
