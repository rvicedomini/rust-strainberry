use std::fmt;

use crate::seq::SeqInterval;

use super::haplotree::Snv;

// pub struct PhasesetRegion {
//     tid: usize,
//     beg: usize,
//     end: usize,
// }

// impl PhasesetRegion {
//     fn length(&self) -> usize { self.end - self.beg }
// }

#[derive(Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct HaplotypeId {
    pub tid: usize,
    pub beg: usize,
    pub end: usize,
    pub hid: usize
}

#[derive(Clone)]
pub struct Haplotype {
    hid: usize,
    tid: usize,
    vars: Vec<Snv>,
    offset: usize,
}

impl Haplotype {

    pub fn new(hid:usize, tid:usize, vars:Vec<Snv>) -> Haplotype {
        Haplotype::with_offset(hid,tid,vars,0)
    }

    pub fn with_offset(hid:usize, tid:usize, vars:Vec<Snv>, offset:usize) -> Haplotype {
        Haplotype { hid, tid, vars, offset }
    }

    pub fn hid(&self) -> usize { self.hid }
    pub fn set_hid(&mut self, hid:usize) { self.hid = hid }

    pub fn tid(&self) -> usize { self.tid }

    pub fn beg(&self) -> usize { self.vars[self.offset].pos }
    pub fn end(&self) -> usize { self.vars.last().unwrap().pos + 1 }

    pub fn uid(&self) -> HaplotypeId {
        HaplotypeId { tid:self.tid(), beg:self.beg(), end:self.end(), hid:self.hid() }
    }

    pub fn size(&self) -> usize { self.vars.len() - self.offset }
    pub fn raw_size(&self) -> usize { self.vars.len() }

    // pub fn get(&self, idx:usize) -> &Snv {
    //     let idx = self.offset + idx;
    //     &self.vars[idx]
    // }

    // pub fn get_mut(&mut self, idx:usize) -> &mut Snv {
    //     let idx = self.offset + idx;
    //     &mut self.vars[idx]
    // }

    pub fn variants(&self) -> &[Snv] { &self.vars[self.offset..] }
    pub fn raw_variants(&self) -> &Vec<Snv> { &self.vars }

    pub fn first(&self) -> &Snv { &self.vars[self.offset] }
    pub fn first_pos(&self) -> usize { self.first().pos }
    // pub fn first_nuc(&self) -> u8 { self.first().nuc }

    pub fn last(&self) -> &Snv { &self.vars[self.vars.len()-1] }
    pub fn last_mut(&mut self) -> &mut Snv { self.vars.last_mut().unwrap() }
    pub fn last_pos(&self) -> usize { self.last().pos }
    pub fn last_nuc(&self) -> u8 { self.last().nuc }

    pub fn region(&self) -> SeqInterval {
        SeqInterval{ tid:self.tid, beg:self.beg(), end:self.end() }
    }

    // pub fn append(&mut self, mut other: Haplotype) {
    //     self.vars.extend(other.vars.drain(other.offset..));
    // }

    pub fn push(&mut self, snv:Snv) {
        self.vars.push(snv);
    }

    pub fn split_off(&mut self, idx:usize, hid:usize, dist:usize) -> Haplotype {
        assert!(self.offset < idx && idx < self.vars.len());
        let target_pos = self.vars[idx].pos;
        let p = self.vars.partition_point(|snv| snv.pos < target_pos.saturating_sub(dist));
        let vars = self.vars[p..].to_vec();
        let offset = idx-p;
        self.vars.truncate(idx);
        Haplotype::with_offset(hid, self.tid, vars, offset)
    }

    // pub fn split_from_position(&mut self, pos:usize, hid:usize, dist:usize) -> Haplotype {
    //     let idx = self.vars.partition_point(|snv| snv.pos < pos);
    //     self.split_off(idx, hid, dist)
    // }

    // pub fn split_until_position(&mut self, pos:usize, hid:usize, dist:usize) -> Haplotype {
    //     let idx = self.vars.partition_point(|snv| snv.pos <= pos);
    //     self.split_off(idx, hid, dist)
    // }

    // pub fn trim(&mut self) {
    //     todo!()
    // }

    pub fn seq_string(&self) -> String {
        let context_str = String::from_iter(self.vars[..self.offset].iter().map(|snv| snv.nuc as char));
        let actual_str = String::from_iter(self.vars[self.offset..].iter().map(|snv| snv.nuc as char));
        format!("{}|{}", context_str, actual_str)
    }

}


impl fmt::Display for Haplotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "h{}: {} ({}:{}..={})", self.hid(), self.seq_string(), self.tid(), self.first_pos(), self.last_pos())
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