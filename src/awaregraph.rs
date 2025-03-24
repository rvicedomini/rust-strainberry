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

    pub fn from_alignments(a:&AwareAlignment, b:&AwareAlignment) -> EdgeKey {
        EdgeKey::new(
            a.aware_id,
            flip_strand(a.strand),
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
fn canonical_edgekey(edge_key:&EdgeKey) -> Cow<'_, EdgeKey> {
    if edge_key.is_canonical() {
        return Cow::Borrowed(edge_key)
    }
    let mut edge_key = edge_key.clone();
    edge_key.flip();
    Cow::Owned(edge_key)
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
        let edge_key: EdgeKey = EdgeKey::from_alignments(a, b);
        let edge = self.get_biedge_or_create(&edge_key);
        edge.observations += 1;
        edge.gaps.push((b.query_beg as i32) - (a.query_end as i32));
        // gapseq = self.read_dict[a.query_id].sequence[a.query_end:b.query_start]
        // edge.gapseq[(key[0],key[1])] = gapseq
        // edge.gapseq[(key[2],key[3])] = reverse_complement(gapseq)
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

    fn get_transitive_or_create(&mut self, edge_key:&EdgeKey) -> &mut BiEdge {
        let edge_key= canonical_edgekey(edge_key).into_owned();
        let edge_id = *self.transitives.entry(edge_key).or_insert_with_key(|edge_key| {
            let node_from = self.nodes.get_mut(&edge_key.id_from).unwrap();
            node_from.transitives.push(edge_key.clone());
            let node_to = self.nodes.get_mut(&edge_key.id_to).unwrap();
            node_to.transitives.push(flip_edgekey(edge_key));
            let new_edge_id = self.edge_data.len();
            self.edge_data.push(BiEdge::new(edge_key.clone()));
            new_edge_id
        });
        unsafe { self.edge_data.get_mut(edge_id).unwrap_unchecked() }
    }

    fn remove_nodes(&mut self, node_ids:&[usize]) {
        for nid in node_ids {
            let biedge_keys: TinyVec<[EdgeKey;10]> = self.nodes
                .get_mut(nid).unwrap()
                .edges.iter().cloned()
                .collect();
            self.remove_biedges_from(&biedge_keys);
            // let tedge_keys_iter = self.nodes[id].transitives.iter();
            // self.remove_tedges_from(biedge_keys_iter);
            self.nodes.remove(nid);
        }
    }

    fn remove_biedge(&mut self, edge_key:&EdgeKey) {
        let edge_key = canonical_edgekey(edge_key);
        self.edges.remove(&edge_key);
    }

    fn remove_biedges_from(&mut self, edge_keys:&[EdgeKey]) {
        for ekey in edge_keys {
            let node = self.nodes.get_mut(&ekey.id_from).unwrap();
            node.edges.retain(|key| key != ekey);
            let ekey = flip_edgekey(ekey);
            let node = self.nodes.get_mut(&ekey.id_from).unwrap();
            node.edges.retain(|key| key != &ekey);
            self.remove_biedge(&ekey);
        }
    }

    fn node_degree(&self, node_id:usize, dir:u8) -> usize {
        self.nodes[&node_id].edges.iter()
            .filter(|&key| key.strand_from == dir)
            .count()
    }

    pub fn remove_weak_edges(&mut self, min_obs:usize) {
        let weak_edges = self.edges.values()
            .filter_map(|&edge_id| {
                let edge = self.get_biedge_idx(edge_id);
                if edge.observations < min_obs {
                    let a_degree = self.node_degree(edge.key.id_from, edge.key.strand_from);
                    let b_degree = self.node_degree(edge.key.id_to, edge.key.strand_to);
                    if a_degree != 1 && b_degree != 1 { // delete only edges that would not disconnect another node
                        return Some(&edge.key)
                    }
                }
                None
            }).cloned().collect_vec();
        self.remove_biedges_from(&weak_edges);
    }

    fn clear_transitive_edges(&mut self) {
        self.nodes.values_mut().for_each(|node| node.transitives.clear());
        self.transitives.clear();
    }

    pub fn add_bridges(&mut self, aware_alignments:&FxHashMap<String,Vec<AwareAlignment>>) -> usize {
        self.clear_transitive_edges();
        // add bridges from aware alignments
        for (_query_name, alignments) in aware_alignments.iter() {
            let mut first_indices: Vec<usize> = vec![];
            let mut last_indices: Vec<usize> = vec![];
            // find bifurcation nodes
            for (idx,(a,b)) in alignments.iter().tuple_windows().enumerate() {
                let EdgeKey { id_from, strand_from, id_to, strand_to } = EdgeKey::from_alignments(a, b);
                if a.aware_id == b.aware_id || !self.nodes.contains_key(&id_from) || !self.nodes.contains_key(&id_to) {
                    continue
                }
                let a_degree = self.node_degree(id_from, strand_from);
                let b_degree = self.node_degree(id_to, strand_to);
                if b_degree > 1 { first_indices.push(idx); }
                if a_degree > 1 { last_indices.push(idx+1); }
            }
            // add transitive edges
            for first_idx in first_indices {
                let first_ctg_id = alignments[first_idx].aware_id;
                let ip = last_indices.partition_point(|&val| val <= first_idx+1);
                for &last_idx in &last_indices[ip..last_indices.len()] {
                    let last_ctg_id = alignments[last_idx].aware_id;
                    if first_ctg_id != last_ctg_id {
                        let a = &alignments[first_idx];
                        let b = &alignments[last_idx];
                        let tedge_key = EdgeKey::from_alignments(a, b);
                        let tedge = self.get_transitive_or_create(&tedge_key);
                        tedge.observations += 1;
                        tedge.gaps.push((b.query_beg as i32) - (a.query_end as i32));
                        // # gapseq = self.read_dict[query_name].sequence[a.query_end:b.query_start]
                        // # tedge.gapseq[(tedge_key[0],tedge_key[1])] = gapseq # (query_name,a.query_end,b.query_start,'+')
                        // # tedge.gapseq[(tedge_key[2],tedge_key[3])] = reverse_complement(gapseq) # (query_name,a.query_end,b.query_start,'-')
                    }
                }
            }
        }
        self.transitives.len()
    }

    pub fn write_gfa(&self, gfa_path:PathBuf, target_names:&[String]) -> std::io::Result<()> {
        let mut gfa = crate::utils::get_file_writer(gfa_path.as_path());
        gfa.write_all(b"H\tVN:Z:1.0\n")?;
        for (node_id, node) in self.nodes.iter() {
            let node_ctg = node.ctg;
            let node_name = format!("{}_{}-{}_h{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), node_ctg.hid().unwrap_or_default(), node_id);
            let node_line = format!("S\t{}\t*\tLN:i:{}\tdp:i:{}\n", node_name, node_ctg.length(), node_ctg.depth() as usize);
            gfa.write_all(node_line.as_bytes())?;
        }
        for edge_id in self.edges.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let node_ctg = self.nodes[&id_from].ctg;
            let name_from = format!("{}_{}-{}_h{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), node_ctg.hid().unwrap_or_default(), id_from);
            let node_ctg = self.nodes[&id_to].ctg;
            let name_to = format!("{}_{}-{}_h{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), node_ctg.hid().unwrap_or_default(), id_to);
            let edge_line = format!("L\t{}\t{}\t{}\t{}\t0M\tRC:i:{}\n", name_from, flip_strand(strand_from) as char, name_to, strand_to as char, edge.observations );
            gfa.write_all(edge_line.as_bytes())?;
        }
        Ok(())
    }

    pub fn write_dot(&self, dot_path:PathBuf) -> std::io::Result<()> {
        let mut dot = crate::utils::get_file_writer(dot_path.as_path());
        dot.write_all(b"digraph \"\" {\n")?;
        dot.write_all(b"\tgraph [rankdir=LR, splines=true];\n")?;
        for (node_id, node) in self.nodes.iter() {
            let node_depth = node.ctg.depth() as usize;
            let fillcolor = ["white","orange"][node.ctg.is_phased() as usize];
            let node_line = format!("\t{node_id}\t[label=\"{node_id} ({node_depth}X)\", style=filled, fillcolor={fillcolor}];\n");
            dot.write_all(node_line.as_bytes())?;
        }
        for edge_id in self.edges.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let arrow_tail = if strand_from == b'+' { "normal" } else { "inv" };
            let arrow_head = if strand_to == b'+' { "normal" } else { "inv" };
            let edge_gap = if edge.gaps.is_empty() { 0 } else { *edge.gaps.iter().sorted_unstable().nth(edge.gaps.len()/2).unwrap() };
            let edge_label = format!("{}/{}bp", edge.observations, edge_gap);
            let edge_line = format!("\t{id_from} -> {id_to}\t[arrowtail={arrow_tail}, arrowhead={arrow_head}, dir=both, label=\"{edge_label}\"];\n");
            dot.write_all(edge_line.as_bytes())?;
        }
        for edge_id in self.transitives.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let arrow_tail = if strand_from == b'+' { "normal" } else { "inv" };
            let arrow_head = if strand_to == b'+' { "normal" } else { "inv" };
            let edge_gap = if edge.gaps.is_empty() { 0 } else { *edge.gaps.iter().sorted_unstable().nth(edge.gaps.len()/2).unwrap() };
            let edge_label = format!("{}/{}bp", edge.observations, edge_gap);
            let edge_line = format!("\t{id_from} -> {id_to}\t[arrowtail={arrow_tail}, arrowhead={arrow_head}, dir=both, color=\"red\", label=\"{edge_label}\"];\n");
            dot.write_all(edge_line.as_bytes())?;
        }
        dot.write_all(b"}\n")
    }

}

