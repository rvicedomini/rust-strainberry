use std::borrow::Cow;
use tinyvec::{tiny_vec,TinyVec};

use crate::awarecontig::{AwareContig, AwareAlignment};

#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct EdgeKey {
    pub id_from: usize,
    pub strand_from: u8,
    pub id_to: usize,
    pub strand_to: u8,
}

impl EdgeKey {

    pub fn new(id_from:usize, strand_from:u8, id_to:usize, strand_to:u8) -> Self {
        Self { id_from, strand_from, id_to, strand_to }
    }

    pub fn from_alignments(a:&AwareAlignment, b:&AwareAlignment) -> EdgeKey {
        EdgeKey::new(
            a.aware_id,
            crate::utils::flip_strand(a.strand),
            b.aware_id,
            b.strand
        )
    }

    pub fn flip(&mut self) {
        std::mem::swap(&mut self.id_from, &mut self.id_to);
        std::mem::swap(&mut self.strand_from, &mut self.strand_to);
    }

    pub fn is_canonical(&self) -> bool {
        self.id_from < self.id_to || (self.id_from == self.id_to && self.strand_from <= self.strand_to)
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
        write!(f, "EdgeKey({},{},{},{})", self.id_from, self.strand_from as char, self.id_to, self.strand_to as char)
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
    pub observations: usize,
    pub gaps: Vec<i32>,
    pub seq_desc: Vec<AwareContig>,
}

impl BiEdge {

    pub fn new(key:EdgeKey) -> Self {
        Self {
            key,
            observations:0,
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
