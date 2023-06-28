use std::mem::MaybeUninit;
// use std::path::Path;

// use itertools::Itertools;

// use rust_htslib::bam;
// use rust_htslib::bam::Read;

// #[derive(Clone)]
// pub struct Region {
//     pub tid: usize,
//     pub chrom: String,
//     pub start_pos: usize,
//     pub end_pos: usize,
// }

// pub fn get_chromosome_regions(bam_path: &Path) -> Vec<Region> {
//     let bam = bam::Reader::from_path(bam_path).unwrap();
//     let header_view = bam.header();
//     let tid2name: Vec<&[u8]> = header_view.target_names();

//     tid2name.iter()
//         .enumerate()
//         .map(|(tid,&name)| Region {
//             tid: tid,
//             chrom: String::from_utf8_lossy(name).to_string(),
//             start_pos: 0,
//             end_pos: header_view.target_len(tid as u32).unwrap() as usize
//         }).collect_vec()
// }

pub fn get_maxrss() -> f64 {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as f64 / (1024.0 * 1024.0)
}