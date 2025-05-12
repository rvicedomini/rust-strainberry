use std::borrow::Cow;
use tinyvec::{tiny_vec,TinyVec};

use crate::awarecontig::{AwareContig, AwareAlignment};

#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct EdgeKey {
    src: usize,
    dst: usize,
}

impl EdgeKey {

    pub fn new(src_id:usize, src_dir:u8, dst_id:usize, dst_dir:u8) -> Self {
        let src = (src_id << 1) | (src_dir != b'+') as usize;
        let dst = (dst_id << 1) | (dst_dir != b'+') as usize;
        Self { src, dst }
    }

    pub fn from_alignments(a:&AwareAlignment, b:&AwareAlignment) -> EdgeKey {
        EdgeKey::new(
            a.aware_id,
            crate::seq::flip_strand(a.strand),
            b.aware_id,
            b.strand
        )
    }

    #[inline(always)]
    pub fn src_id(&self) -> usize {
        self.src >> 1
    }

    #[inline(always)]
    pub fn src_dir(&self) -> u8 {
        if (self.src & 1) == 0 { b'+' } else { b'-' }
    }

    #[inline(always)]
    pub fn src(&self) -> (usize,u8) {
        (self.src_id(), self.src_dir())
    }

    #[inline(always)]
    pub fn dst_id(&self) -> usize {
        self.dst >> 1
    }

    #[inline(always)]
    pub fn dst_dir(&self) -> u8 {
        if (self.dst & 1) == 0 { b'+' } else { b'-' }
    }

    #[inline(always)]
    pub fn dst(&self) -> (usize,u8) {
        (self.dst_id(), self.dst_dir())
    }

    #[inline(always)]
    pub fn unpack(&self) -> (usize,u8,usize,u8) {
        (self.src_id(), self.src_dir(), self.dst_id(), self.dst_dir())
    }

    #[inline(always)]
    pub fn flip(&mut self) {
        std::mem::swap(&mut self.src, &mut self.dst);
    }

    #[inline(always)]
    pub fn is_canonical(&self) -> bool {
        self.src <= self.dst
    }

    // returns whether the key was already canonical
    pub fn canonicalize(&mut self) -> bool {
        if !self.is_canonical() {
            self.flip();
            return false
        }
        true
    }

}

impl std::fmt::Display for EdgeKey {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "EdgeKey({},{},{},{})", self.src_id(), self.src_dir() as char, self.dst_id(), self.dst_dir() as char)
    }
}


#[inline(always)]
pub fn flip_edgekey(key:&EdgeKey) -> EdgeKey {
    let mut key = *key;
    key.flip();
    key
}


#[inline(always)]
pub fn canonical_edgekey(edge_key:&EdgeKey) -> Cow<'_, EdgeKey> {
    if edge_key.is_canonical() {
        return Cow::Borrowed(edge_key)
    }
    let mut edge_key = *edge_key;
    edge_key.flip();
    Cow::Owned(edge_key)
}

#[derive(Debug)]
pub struct BiEdge {
    pub key: EdgeKey, // TODO: check whether I really need to store the key in the struct
    pub nb_reads: usize,
    pub min_shared_snvs: usize,
    pub gaps: Vec<i32>,
    pub seq_desc: Vec<AwareContig>,
}

impl BiEdge {

    pub fn new(key:EdgeKey) -> Self {
        Self {
            key,
            nb_reads:0,
            min_shared_snvs:0,
            gaps: Vec::new(),
            seq_desc: Vec::new()
        }
    }
}


#[derive(Debug)]
pub struct Node {
    pub id: usize,
    pub ctg: AwareContig,
    pub edges: TinyVec<[EdgeKey;10]>,
    pub transitives: TinyVec<[EdgeKey;10]>,
}

impl Node {

    pub fn new(id:usize, ctg: AwareContig) -> Self {
        Self {
            id,
            ctg,
            edges: tiny_vec![],
            transitives: tiny_vec![]
        }
    }
}
