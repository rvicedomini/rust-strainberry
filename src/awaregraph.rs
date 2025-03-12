use rustc_hash::FxHashMap;
use tinyvec::{tiny_vec,TinyVec};

use crate::awarecontig::AwareContig;


#[inline(always)]
fn flip_strand(strand:u8) -> u8 {
    assert!(strand == b'+' || strand == b'-');
    // b'+': 00101|01|1
    // b'-': 00101|10|1
    strand ^ 6
}


#[derive(Debug)]
struct Node<'a> {
    id: usize,
    ctg: &'a AwareContig,
    edges: TinyVec<[EdgeKey;10]>,
    transitives: TinyVec<[EdgeKey;10]>,
}

impl<'a> Node<'a> {

    pub fn new(id:usize, ctg:&'a AwareContig) -> Self {
        Self { id, ctg, edges: tiny_vec![], transitives: tiny_vec![] }
    }
}


// TODO: move BiEdge/TransEdge in a sub-module?

type EdgeKey = (usize, u8, usize, u8); // (id_from, strand_from, id_to, strand_to)

#[inline(always)]
fn flip_edgekey(key:EdgeKey) -> EdgeKey {
    (key.2, key.3, key.0, key.1)
}

#[inline(always)]
fn canon_edgekey(key:EdgeKey) -> EdgeKey {
    if (key.0,key.1) <= (key.2,key.3) {
        key
    } else {
        flip_edgekey(key)
    }
}

#[derive(Debug)]
struct BiEdge {
    key: EdgeKey, 
    observations: usize,
    gaps: Vec<i32>,
    // self.gapseq = {} # defaultdict(list)
}

impl BiEdge {

    pub fn new(key:EdgeKey) -> Self {
        Self {
            key,
            observations:0,
            gaps:vec![]
        }
    }
}


#[derive(Debug, Default)]
pub struct AwareGraph<'a> {
    next_node_id: usize,
    nodes: FxHashMap<usize,Node<'a>>,
    edges: FxHashMap<EdgeKey,BiEdge>,
    transitives: FxHashMap<EdgeKey,BiEdge>,
}


