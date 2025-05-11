use std::path::Path;

use ahash::AHashMap as HashMap;
use needletail::Sequence;

// from ffforf: https://github.com/jguhlin/ffforf/blob/master/src/lib.rs
#[inline(always)]
fn complement(nuc: u8) -> u8 {
    if nuc != b'N' {
        if nuc & 2 != 0 {
            nuc ^ 4
        } else {
            nuc ^ 21
        }
    } else {
        nuc
    }
}


pub fn revcomp_inplace(seq: &mut [u8]) {
    seq.reverse();
    seq.iter_mut().for_each(|nuc| { *nuc = complement(*nuc) });
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut rev = seq.to_vec();
    revcomp_inplace(&mut rev);
    rev
}


pub fn load_sequences(fastx_path: &Path) -> (Vec<Vec<u8>>,HashMap<String,usize>) {
    let mut sequences = Vec::new();
    let mut name_index = HashMap::new();

    let mut reader = needletail::parse_fastx_file(fastx_path).expect("Cannot open file");
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let seq = record.normalize(false).to_vec();
        let name = record.id().split(|b| b.is_ascii_whitespace()).next().unwrap();
        let name = String::from_utf8_lossy(name).to_string();
        name_index.insert(name, sequences.len());
        sequences.push(seq);
    }

    (sequences, name_index)
}
