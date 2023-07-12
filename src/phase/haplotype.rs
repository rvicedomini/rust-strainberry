use itertools::Itertools;

pub struct PhasesetRegion {
    tid: usize,
    beg: usize,
    end: usize,
}

impl PhasesetRegion {
    fn length(&self) -> usize { self.end - self.beg }
}

pub struct Haplotype {
    hid: usize,
    tid: usize,
    positions: Vec<usize>,
    nucleotides: Vec<u8>,
    offset: usize,
}

impl Haplotype {

    pub fn hid(&self) -> usize { self.hid }

    pub fn tid(&self) -> usize { self.tid }

    pub fn beg(&self) -> usize { self.positions[self.offset] }

    pub fn end(&self) -> usize { self.positions.last().unwrap() + 1 }

    pub fn size(&self) -> usize { self.nucleotides.len() - self.offset }

    pub fn at(&self, idx:usize) -> (usize,u8) {
        let idx = self.offset + idx;
        (self.positions[idx], self.nucleotides[idx])
    }

    pub fn phaseset_region(&self) -> PhasesetRegion {
        PhasesetRegion{ tid:self.tid, beg:self.beg(), end:self.end() }
    }

    pub fn extend(&mut self, other:Haplotype) {
        self.positions.extend(&other.positions[other.offset..]);
        self.nucleotides.extend(&other.nucleotides[other.offset..])
    }

    pub fn append(&mut self, pos:usize, nuc:u8) {
        self.positions.push(pos);
        self.nucleotides.push(nuc);
    }

    pub fn split_off(&mut self, idx:usize, hid:usize, dist:usize) -> Haplotype {
        assert!(self.offset < idx && idx < self.positions.len());
        let tid = self.tid;
        let p_idx = self.positions.partition_point(|&pos| pos < self.positions[idx]-dist);
        let positions = self.positions.split_off(idx);
        let nucleotides = self.nucleotides.split_off(idx);
        Haplotype { hid, tid, positions, nucleotides, offset: idx - p_idx }
    }

    pub fn split_from_position(&mut self, pos:usize, hid:usize, dist:usize) -> Haplotype {
        let idx = self.positions.partition_point(|&p| p < pos);
        self.split_off(idx, hid, dist)
    }

    pub fn split_until_position(&mut self, pos:usize, hid:usize, dist:usize) -> Haplotype {
        let idx = self.positions.partition_point(|&p| p <= pos);
        self.split_off(idx, hid, dist)
    }

    pub fn trim(&mut self) {
        todo!()
    }

}

// def trim(self, variant_pos, lookback):
//     # left_trim
//     it = (i for i in range(len(self.positions)) if self.positions[i] in variant_pos)
//     first_i = next(it,None)
//     if first_i is None: # no variant positions -> delete haplotype
//         self.positions = []
//         self.sequence = ''
//         self.start_idx = 0
//         return self
//     # trim positions that are too far from the first variant position
//     var_pos = self.positions[first_i]
//     cut_i = first_i-1
//     while cut_i >= 0 and var_pos-self.positions[cut_i]+1 <= lookback:
//         cut_i -= 1
//     if cut_i >= 0:
//         self.positions = self.positions[cut_i+1:]
//         self.sequence = self.sequence[cut_i+1:]
//         self.start_idx = max(0,self.start_idx-(cut_i+1))
//     # right trim
//     it = (i for i in reversed(range(self.size())) if self.positions[self.start_idx+i] in variant_pos)
//     last_i = next(it,None)
//     if last_i is None: # no variant positions -> delete haplotype
//         self.positions = []
//         self.sequence = ''
//         self.start_idx = 0
//         return self
//     # trim positions that are too far from the last variant position
//     var_pos = self.positions[self.start_idx+last_i]
//     cut_i = self.start_idx + last_i + 1
//     while cut_i < len(self.positions) and self.positions[cut_i]-var_pos+1 <= lookback:
//         cut_i += 1
//     if cut_i < len(self.positions):
//         self.positions = self.positions[:cut_i]
//         self.sequence = self.sequence[:cut_i]
//     return self