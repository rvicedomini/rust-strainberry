use std::io::BufRead;
use std::path::Path;

use anyhow::{bail, Context, Result};
use ahash::AHashMap as HashMap;
use itertools::Itertools;
use rust_htslib::bam::record::{Cigar,CigarString};


const ERROR_THRESHOLD: f64 = 0.3;


#[derive(Debug, Clone)]
pub struct Alignment {
    pub query_idx: usize,
    pub query_len: usize,
    pub query_beg: usize,
    pub query_end: usize,
    pub strand: u8,
    pub target_idx: usize,
    pub target_len: usize,
    pub target_beg: usize,
    pub target_end: usize,
    pub matches: usize, // exact matches
    pub mapping_length: usize,
    pub mapq: u8,
    pub length: usize,
    pub identity: f64,
    pub cigar: CigarString,
    pub breaking_points: Vec<(usize,usize)>
}

impl Alignment {

    fn parse_paf_line(line: &str, ref_index: &HashMap<String,usize>, read_index: &HashMap<String,usize>) -> Result<Self> {
        
        // let cols: [&str; 13] = line.splitn(13, '\t').collect_array().context("cannot parse PAF line")?;
        let cols = line.split('\t').collect_vec();
        if cols.len() < 12 { bail!("cannot parse PAF line (mission fields)") }
        
        let query_idx = *read_index.get(cols[0]).with_context(|| format!("missing read with id: {}", cols[0]))?;
        let [query_len, query_beg, query_end] = cols[1..4].iter().map(|s|s.parse::<usize>().unwrap()).collect_array().unwrap();
        
        let strand = match cols[4] {
            "+" => b'+',
            "-" => b'-',
            _ => bail!("unrecognised strand field: {}", cols[4])
        };
        
        let target_idx = ref_index[cols[5]];
        let [target_len, target_beg, target_end] = cols[6..9].iter().map(|s|s.parse::<usize>().unwrap()).collect_array().unwrap();
        
        let matches = cols[9].parse().unwrap();
        let mapping_length = cols[10].parse().unwrap();
        let mapq = cols[11].parse().unwrap();
        let identity = 100.0 * (matches as f64) / (mapping_length as f64);

        let length = std::cmp::max(query_end-query_beg, target_end-target_beg);
        
        Ok(Alignment {
            query_idx,
            query_len,
            query_beg,
            query_end,
            strand,
            target_idx,
            target_len,
            target_beg,
            target_end,
            matches,
            mapping_length,
            mapq,
            length,
            identity,
            cigar: CigarString(Vec::new()),
            breaking_points: Vec::new(),
        })
    }


    pub fn find_breaking_points(&mut self, ref_sequences: &[Vec<u8>], read_sequences: &[Vec<u8>], window_len: usize) -> Result<()> {

        if !self.breaking_points.is_empty() {
            return Ok(())
        }

        if self.cigar.is_empty() { // compute alignment

            use edlib_rs::edlibrs;

            let tseq = &ref_sequences[self.target_idx][self.target_beg..self.target_end];
            let mut qseq = read_sequences[self.query_idx][self.query_beg..self.query_end].to_vec();
            if self.strand == b'-' { crate::racon::sequence::revcomp_inplace(&mut qseq); }

            let ed_cfg = edlibrs::EdlibAlignConfigRs::new(
                -1,
                edlibrs::EdlibAlignModeRs::EDLIB_MODE_NW,
                edlibrs::EdlibAlignTaskRs::EDLIB_TASK_PATH,
                &[]
            );

            let align_res = edlibrs::edlibAlignRs(&qseq, &tseq, &ed_cfg);
            if align_res.status != edlibrs::EDLIB_STATUS_OK {
                bail!("edlib: unable to align query {} against target {}", self.query_idx, self.target_idx);
            }

            let alignment = align_res.alignment.as_ref().unwrap();
            let cigar = edlibrs::edlibAlignmentToCigarRs(alignment, &edlibrs::EdlibCigarFormatRs::EDLIB_CIGAR_STANDARD);
            self.cigar = parse_cigar_bytes(cigar.as_bytes());
        }
    
        self.find_breaking_points_from_cigar(window_len);

        Ok(())
    }


    pub fn find_breaking_points_from_cigar(&mut self, window_len: usize) {

        // find breaking points from cigar
        let window_ends: Vec<usize> = (0..self.target_end)
            .step_by(window_len)
            .filter_map(|i| if i > self.target_beg { Some(i-1) } else { None })
            .chain(std::iter::once(self.target_end-1))
            .collect();

        let mut w = 0;
        let mut found_first_match = false;
        let mut first_match: (usize,usize) = (0,0);
        let mut last_match: (usize,usize) = (0,0);

        let mut q_ptr: i32 = if self.strand == b'+' { self.query_beg } else { self.query_len-self.query_end } as i32 - 1;
        let mut t_ptr: i32 = self.target_beg as i32 - 1;

        for cigar_op in &self.cigar {
            match *cigar_op {
                Cigar::Match(num_bases) | Cigar::Equal(num_bases) | Cigar::Diff(num_bases) => {
                    for _ in 0..num_bases {
                        q_ptr += 1; t_ptr += 1;
                        if !found_first_match {
                            found_first_match = true;
                            first_match = (t_ptr as usize, q_ptr as usize);
                        }
                        last_match = (t_ptr as usize + 1, q_ptr as usize + 1);
                        if t_ptr as usize == window_ends[w] {
                            if found_first_match {
                                self.breaking_points.push(first_match);
                                self.breaking_points.push(last_match);
                            }
                            found_first_match = false;
                            w += 1;
                        }
                    }
                },
                Cigar::Ins(num_bases) | Cigar::SoftClip(num_bases) => {
                    q_ptr += num_bases as i32;
                }
                Cigar::Del(num_bases) | Cigar::RefSkip(num_bases) => {
                    for _ in 0..num_bases {
                        t_ptr += 1;
                        if t_ptr as usize == window_ends[w] {
                            if found_first_match {
                                self.breaking_points.push(first_match);
                                self.breaking_points.push(last_match);
                            }
                            found_first_match = false;
                            w += 1;
                        }
                    }
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {
                    continue;
                }
            }
        }
    }
}


pub fn load_paf_alignments(paf_path: &Path, ref_index: &HashMap<String,usize>, read_index: &HashMap<String,usize>) -> Result<Vec<Alignment>> {

    let mut alignments = Vec::new();
    let reader = crate::utils::get_file_reader(paf_path);
    for line in reader.lines().map_while(Result::ok) {
        let line = line.trim();
        if line.is_empty() { continue }
        let paf = Alignment::parse_paf_line(line, ref_index, read_index)?;
        alignments.push(paf);
    }

    Ok(alignments)
}


pub fn filter_alignments(alignments: Vec<Alignment>) -> Vec<Alignment> {

    let mut retained: HashMap<usize, Alignment> = HashMap::new();

    for a in alignments {
        let a_err = 1.0 - std::cmp::min(a.query_end-a.query_beg, a.target_end-a.target_beg) as f64 / a.length  as f64;
        if a_err > ERROR_THRESHOLD {
            continue;
        }
        
        retained.entry(a.query_idx)
            .and_modify(|e| {
                if a.length > e.length { *e = a.clone(); }
            }).or_insert(a);
    }

    retained.into_values().collect()
}


pub fn parse_cigar_bytes(cigar: &[u8]) -> CigarString {
    CigarString::try_from(cigar)
        .expect("Unable to parse cigar string.")
}

