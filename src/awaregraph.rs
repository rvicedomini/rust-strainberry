use std::borrow::Cow;
use std::path::PathBuf;

use itertools::Itertools;
use rustc_hash::FxHashMap;
use tinyvec::{tiny_vec,TinyVec};

use crate::awarecontig::{AwareContig,AwareAlignment};


#[inline(always)]
fn flip_strand(strand:u8) -> u8 {
    assert!(strand == b'+' || strand == b'-');
    // b'+': 00101|01|1
    // b'-': 00101|10|1
    //    6: 00000|11|0
    strand ^ 6
}


#[derive(Debug)]
struct Node<'a> {
    pub id: usize,
    pub ctg: &'a AwareContig,
    pub edges: TinyVec<[EdgeKey;10]>,
    pub transitives: TinyVec<[EdgeKey;10]>,
}

impl<'a> Node<'a> {

    pub fn new(id:usize, ctg:&'a AwareContig) -> Self {
        Self { id, ctg, edges: tiny_vec![], transitives: tiny_vec![] }
    }
}


// TODO: move BiEdge in a sub-module?

//type EdgeKey = (usize, u8, usize, u8); // (id_from, strand_from, id_to, strand_to)

#[derive(Debug, Default, Clone, Hash, PartialEq, Eq)]
struct EdgeKey {
    pub id_from: usize,
    pub strand_from: u8,
    pub id_to: usize,
    pub strand_to: u8,
}

impl EdgeKey {

    pub fn new(id_from:usize, strand_from:u8, id_to:usize, strand_to:u8) -> Self {
        Self { id_from, strand_from, id_to, strand_to }
    }

    pub fn flip(&mut self) {
        std::mem::swap(&mut self.id_from, &mut self.id_to);
        std::mem::swap(&mut self.strand_from, &mut self.strand_to);
    }

    pub fn is_canonical(&self) -> bool {
        self.id_from < self.id_to || (self.id_from == self.id_to && self.strand_from <= self.strand_to)
    }

    pub fn canonicalize(&mut self) {
        if !self.is_canonical() {
            self.flip();
        }
    }

}

#[inline(always)]
fn flip_edgekey(key:&EdgeKey) -> EdgeKey {
    let mut key = key.clone();
    key.flip();
    key
}

#[inline(always)]
fn canonical_edgekey(key:&EdgeKey) -> Cow<'_, EdgeKey> {
    if key.is_canonical() {
        return Cow::Borrowed(key)
    }
    let mut key = key.clone();
    key.flip();
    Cow::Owned(key)
}

// TODO: check whether I really need to store the key in the struct
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
    nodes: FxHashMap<usize,Node<'a>>,
    edges: FxHashMap<EdgeKey,usize>,
    transitives:FxHashMap<EdgeKey,usize>,
    edge_data: Vec<BiEdge>,
    // next_node_id: usize,
}

impl<'a> AwareGraph<'a> {

    pub fn build(aware_contigs:&'a [AwareContig]) -> Self {
        let contigs_iter = aware_contigs.iter().enumerate()
            .map(|(node_id,aware_contig)| (node_id, Node::new(node_id, aware_contig)));
        let nodes = FxHashMap::from_iter(contigs_iter);
        Self {
            nodes,
            edges: FxHashMap::default(),
            transitives: FxHashMap::default(),
            edge_data: vec![]
        }
    }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn nb_edges(&self) -> usize { self.edges.len() }

    pub fn len(&self) -> usize { self.nb_nodes() }
    pub fn is_empty(&self) -> bool { self.nb_nodes() == 0 }

    pub fn add_edges_from_aware_alignments(&mut self, aware_alignments:&FxHashMap<String,Vec<AwareAlignment>>) {
        for alignments in aware_alignments.values() {
            for (a, b) in alignments.iter().tuple_windows() {
                self.add_edge_from_consecutive_alignments(a,b);
            }
        }
    }

    fn add_edge_from_consecutive_alignments(&mut self, a:&AwareAlignment, b:&AwareAlignment) {
        let edge_key: EdgeKey = Self::get_biedge_key(a, b);
        let edge = self.get_biedge_or_create(&edge_key);
        edge.observations += 1;
        edge.gaps.push((b.query_beg as i32) - (a.query_end as i32));
        // gapseq = self.read_dict[a.query_id].sequence[a.query_end:b.query_start]
        // edge.gapseq[(key[0],key[1])] = gapseq
        // edge.gapseq[(key[2],key[3])] = reverse_complement(gapseq)
    }

    fn get_biedge_key(a:&AwareAlignment, b:&AwareAlignment) -> EdgeKey {
        EdgeKey::new(
            a.aware_id,
            flip_strand(a.strand),
            b.aware_id,
            b.strand
        )
    }

    fn contains_edge(&self, edge_key: &EdgeKey) -> bool {
        let edge_key = canonical_edgekey(edge_key);
        self.edges.contains_key(&edge_key)
    }

    fn get_biedge_idx(&self, idx: usize) -> &BiEdge { &self.edge_data[idx] }
    fn get_biedge_idx_mut(&mut self, idx: usize) -> &mut BiEdge { &mut self.edge_data[idx] }

    fn get_biedge(&self, edge_key:&EdgeKey) -> Option<&BiEdge> {
        let edge_key = canonical_edgekey(edge_key);
        let edge_idx = *self.edges.get(&edge_key)?;
        self.edge_data.get(edge_idx)
    }

    fn get_biedge_mut(&mut self, edge_key:&EdgeKey) -> Option<&mut BiEdge> {
        let edge_key = canonical_edgekey(edge_key);
        let edge_idx = *self.edges.get(&edge_key)?;
        self.edge_data.get_mut(edge_idx)
    }

    fn get_biedge_or_create(&mut self, edge_key:&EdgeKey) -> &mut BiEdge {
        let edge_key= canonical_edgekey(edge_key).into_owned();
        let edge_id = *self.edges.entry(edge_key).or_insert_with_key(|edge_key| {
            let node_from = self.nodes.get_mut(&edge_key.id_from).unwrap();
            node_from.edges.push(edge_key.clone());
            let node_to = self.nodes.get_mut(&edge_key.id_to).unwrap();
            node_to.edges.push(flip_edgekey(edge_key));
            let new_edge_id = self.edge_data.len();
            self.edge_data.push(BiEdge::new(edge_key.clone()));
            new_edge_id
        });
        unsafe { self.edge_data.get_mut(edge_id).unwrap_unchecked() }
    }

    pub fn write_gfa(&self, gfa_path:PathBuf) -> std::io::Result<()> {
        let mut gfa = crate::utils::get_file_writer(gfa_path.as_path());
        gfa.write_all(b"H\tVN:Z:1.0\n")?;
        for (node_id, node) in self.nodes.iter() {
            let node_ctg = node.ctg;
            let node_name = format!("{}_{}-{}_h{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), node_ctg.hid().unwrap_or_default(), node_id);
            let node_line = format!("S\t{}\t*\tLN:i:{}\tdp:i:{}\n", node_name, node_ctg.length(), node_ctg.depth() as usize);
            gfa.write_all(node_line.as_bytes())?;
        }
        for edge_id in self.edges.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let node_ctg = self.nodes[&id_from].ctg;
            let name_from = format!("{}_{}-{}_h{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), node_ctg.hid().unwrap_or_default(), id_from);
            let node_ctg = self.nodes[&id_to].ctg;
            let name_to = format!("{}_{}-{}_h{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), node_ctg.hid().unwrap_or_default(), id_to);
            let edge_line = format!("L\t{}\t{}\t{}\t{}\t0M\tRC:i:{}\n", name_from, flip_strand(strand_from) as char, name_to, strand_to as char, edge.observations );
            gfa.write_all(edge_line.as_bytes())?;
        }
        Ok(())
    }


}

