use std::fmt;
use rustc_hash::FxHashSet;

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

#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct HaplotypeId {
    pub tid: usize,
    pub beg: usize,
    pub end: usize,
    pub hid: usize
}

impl fmt::Display for HaplotypeId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}_{}-{}_h{}", self.tid, self.beg, self.end, self.hid)
    }
}


#[derive(Debug, Clone)]
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

    pub fn extend(&mut self, mut other: Haplotype) {
        self.vars.extend(other.vars.drain(other.offset..));
    }

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

    pub fn trim(&mut self, variant_positions:&FxHashSet<(usize,usize)>, lookback:usize) {
        // left trim
        let first_i = (0..self.vars.len()).find(|&i| variant_positions.contains(&(self.tid(),self.vars[i].pos)));
        if first_i.is_none() { // no variant positions -> delete haplotype
            self.vars.clear();
            self.offset = 0;
            return;
        }
        // trim positions that are too far from the first variant position
        let first_i = unsafe { first_i.unwrap_unchecked() }; // if it was None I would have returned in the if-statement above
        let var_pos = self.vars[first_i].pos;
        let cut_i = (0..first_i).rev().find(|i| var_pos - self.vars[*i].pos > lookback);
        if let Some(cut_i) = cut_i {
            self.vars.drain(..=cut_i);
            self.offset = std::cmp::max(0, self.offset-(cut_i+1))
        }
        // right trim
        let last_i = (0..self.size()).rev().find(|&i| variant_positions.contains(&(self.tid(),self.vars[self.offset+i].pos)));
        if last_i.is_none() { // no variant positions -> delete haplotype
            self.vars.clear();
            self.offset = 0;
            return;
        }
        // trim positions that are too far from the last variant position
        let last_i = unsafe { last_i.unwrap_unchecked() }; // if it was None I would have returned in the if-statement above
        let var_pos = self.vars[self.offset+last_i].pos;
        let cut_i = (self.offset+last_i+1..self.vars.len()).find(|i| self.vars[*i].pos - var_pos > lookback);
        if let Some(cut_i) = cut_i {
            self.vars.drain(..cut_i);
        }
    }

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