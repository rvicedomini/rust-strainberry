use anyhow::{Error, Result};

const B2N: &[u8;4] = b"ACGT";

#[derive(Debug, Default, Clone)]
pub struct BitSeq {
    data: Vec<u8>,
    len: usize,
}

impl BitSeq {

    #[inline(always)]
    pub fn from_utf8(seq: &[u8]) -> Self {
        let len = seq.len();
        let nb_bytes = (len/4) + (!len.is_multiple_of(4)) as usize;
        let mut data = vec![0u8;nb_bytes];

        for (i,bp) in seq.iter().enumerate() {
            let base = 0b11 & ((bp >> 2) ^ (bp >> 1)); // convert nucleotide ascii chars "ACGT" to integers (0,1,2,3)
            let shift = 6 - ((i & 0b11) << 1);
            data[i/4] |= base << shift;
        }

        Self { data, len }
    }

    #[inline(always)]
    pub fn len(&self) -> usize { self.len }

    #[inline(always)]
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    #[inline(always)]
    pub fn get(&self, index:usize) -> u8 {
        let shift = 6 - ((index & 0b11) << 1);
        let base = (self.data[index/4] >> shift) & 0b11;
        B2N[base as usize]
    }

    #[inline(always)]
    pub fn as_bytes(&self) -> Vec<u8> { self.subseq(0, self.len()) }

    #[inline(always)]
    pub fn subseq(&self, beg:usize, end:usize) -> Vec<u8> {
        assert!(beg <= end && end <= self.len);
        let mut res = Vec::with_capacity(end-beg);
        (beg..end).for_each(|index| res.push(self.get(index)));
        res
    }

    // #[inline(always)]
    // pub fn revcomp(mut b:u8) -> u8 {
    //     b = ((b & 0xF0) >> 4) | ((b & 0x0F) << 4);
    //     ((b & 0xCC) >> 2) | ((b & 0x33) << 2)
    // }
    
}

impl std::str::FromStr for BitSeq {
    type Err = Error;
    
    #[inline(always)]
    fn from_str(seq: &str) -> Result<Self,Self::Err> {
        Ok(Self::from_utf8(seq.as_bytes()))
    }
}

impl std::fmt::Display for BitSeq {

    #[inline(always)]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let bytes = self.as_bytes();
        write!(f, "{}", std::str::from_utf8(&bytes).unwrap())
    }
}
