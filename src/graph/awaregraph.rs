use std::borrow::Borrow;
use std::path::{Path,PathBuf};

use anyhow::{Context,Result};
use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use itertools::Itertools;
use tinyvec::{TinyVec, ArrayVec};

use crate::awarecontig::{AwareAlignment, AwareContig, SeqType};
use crate::cli::Options;
use crate::graph::asmgraph::{AsmGraph, Link};
use crate::polish::PolishMode;
use crate::seq::{SeqDatabase, SeqInterval};
use crate::bitseq::BitSeq;

use super::biedge::{self, BiEdge, EdgeKey, Node};
use super::junction::Junction;

#[derive(Debug)]
pub struct AwarePath {
    pub nodes: Vec<usize>,
    pub edges: Vec<EdgeKey>,
}


#[derive(Debug, Default)]
pub struct AwareGraph {
    nodes: HashMap<usize,Node>,
    edges: HashMap<EdgeKey,usize>,
    transitives: HashMap<EdgeKey,usize>,
    edge_data: Vec<BiEdge>,
    next_node_id: usize,
}

impl AwareGraph {

    // Prerequisite: aware_contigs must be sorted by interval for the following to work
    pub fn build(aware_contigs:&[AwareContig]) -> Self {

        let contigs_iter = aware_contigs.iter().enumerate()
            .map(|(node_id,aware_contig)| (node_id, Node::new(node_id, *aware_contig)));
        let nodes = HashMap::from_iter(contigs_iter);
        let next_node_id = nodes.len();

        let mut aware_graph = Self {
            nodes,
            edges: HashMap::new(),
            transitives: HashMap::new(),
            edge_data: vec![],
            next_node_id
        };

        // insert edges between contigs adjacent on reference
        aware_contigs.iter().tuple_windows().enumerate()
            .filter(|(_,(a,b))| a.tid() == b.tid() && !a.is_phased() && !b.is_phased())
            .for_each(|(i,(a,b))| {
                let edge_key = EdgeKey::new(i, b'-', i+1, b'+');
                let edge = aware_graph.get_biedge_or_create(&edge_key);
                let interval = SeqInterval{ tid: a.tid(), beg: a.end(), end: b.beg() };
                // strand of edge sequence is defined according to the "direction" of the canonical edge
                assert!(edge_key.is_canonical());
                let aware_contig = AwareContig::new(SeqType::Unphased, interval, b'+', 0);
                edge.seq_desc.push(aware_contig);
            });

        aware_graph
    }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn nb_edges(&self) -> usize { self.edges.len() }

    pub fn len(&self) -> usize { self.nb_nodes() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn add_node(&mut self, ctg: AwareContig) -> usize {
        let node_id = self.next_node_id;
        self.nodes.insert(node_id, Node::new(node_id, ctg));
        self.next_node_id += 1;
        node_id
    }

    pub fn add_edges_from_aware_alignments(&mut self, aware_alignments:&HashMap<usize,Vec<AwareAlignment>>) {
        
        for alignments in aware_alignments.values() {
            let alignments_iter = alignments.iter()
                .tuple_windows()
                .filter(|(a,b)| !a.is_ambiguous && !b.is_ambiguous);
            for (a, b) in alignments_iter {
                self.add_edge_from_consecutive_alignments(a,b);
            }
        }
    }

    fn add_edge_from_consecutive_alignments(&mut self, a:&AwareAlignment, b:&AwareAlignment) {
        let mut edge_key: EdgeKey = EdgeKey::from_alignments(a, b);
        let was_canonical = edge_key.canonicalize();
        let edge = self.get_biedge_or_create(&edge_key);
        edge.nb_reads += 1;
        edge.min_shared_snvs = edge.min_shared_snvs.max(a.nb_shared_snvs.min(b.nb_shared_snvs));
        edge.gaps.push((b.query_beg as i32) - (a.query_end as i32));
        if edge.seq_desc.is_empty() && b.query_beg > a.query_end {
            let interval = SeqInterval{ tid: a.query_idx, beg: a.query_end, end: b.query_beg };
            let strand = if was_canonical { b'+' } else { b'-' };
            let aware_contig = AwareContig::new(SeqType::Read, interval, strand, 0);
            edge.seq_desc.push(aware_contig);
        }
    }

    fn contains_edge(&self, edge_key: &EdgeKey) -> bool {
        let edge_key = biedge::canonical_edgekey(edge_key);
        self.edges.contains_key(&edge_key)
    }

    fn get_biedge_idx(&self, idx: usize) -> &BiEdge { &self.edge_data[idx] }
    fn get_biedge_idx_mut(&mut self, idx: usize) -> &mut BiEdge { &mut self.edge_data[idx] }

    fn get_biedge(&self, edge_key:&EdgeKey) -> Option<&BiEdge> {
        let edge_key = biedge::canonical_edgekey(edge_key);
        let edge_idx = *self.edges.get(&edge_key)?;
        self.edge_data.get(edge_idx)
    }

    fn get_biedge_mut(&mut self, edge_key:&EdgeKey) -> Option<&mut BiEdge> {
        let edge_key = biedge::canonical_edgekey(edge_key);
        let edge_idx = *self.edges.get(&edge_key)?;
        self.edge_data.get_mut(edge_idx)
    }

    fn get_biedge_or_create(&mut self, edge_key:&EdgeKey) -> &mut BiEdge {
        let edge_key= biedge::canonical_edgekey(edge_key).into_owned();
        let edge_id = *self.edges.entry(edge_key).or_insert_with_key(|edge_key| {
            let node_from = self.nodes.get_mut(&edge_key.id_from).unwrap();
            node_from.edges.push(*edge_key);
            let node_to = self.nodes.get_mut(&edge_key.id_to).unwrap();
            node_to.edges.push(biedge::flip_edgekey(edge_key));
            let new_edge_id = self.edge_data.len();
            self.edge_data.push(BiEdge::new(*edge_key));
            new_edge_id
        });
        unsafe { self.edge_data.get_mut(edge_id).unwrap_unchecked() }
    }

    fn get_transitive(&self, edge_key:&EdgeKey) -> Option<&BiEdge> {
        let edge_key = biedge::canonical_edgekey(edge_key);
        let edge_idx = *self.transitives.get(&edge_key)?;
        self.edge_data.get(edge_idx)
    }

    fn get_transitive_mut(&mut self, edge_key:&EdgeKey) -> Option<&mut BiEdge> {
        let edge_key = biedge::canonical_edgekey(edge_key);
        let edge_idx = *self.transitives.get(&edge_key)?;
        self.edge_data.get_mut(edge_idx)
    }

    fn get_transitive_or_create(&mut self, edge_key:&EdgeKey) -> &mut BiEdge {
        let edge_key= biedge::canonical_edgekey(edge_key).into_owned();
        let edge_id = *self.transitives.entry(edge_key).or_insert_with_key(|edge_key| {
            let node_from = self.nodes.get_mut(&edge_key.id_from).unwrap();
            node_from.transitives.push(*edge_key);
            let node_to = self.nodes.get_mut(&edge_key.id_to).unwrap();
            node_to.transitives.push(biedge::flip_edgekey(edge_key));
            let new_edge_id = self.edge_data.len();
            self.edge_data.push(BiEdge::new(*edge_key));
            new_edge_id
        });
        unsafe { self.edge_data.get_mut(edge_id).unwrap_unchecked() }
    }

    fn contains_transitive(&self, edge_key:&EdgeKey) -> bool {
        let edge_key = biedge::canonical_edgekey(edge_key);
        self.transitives.contains_key(&edge_key)
    }

    fn remove_nodes_from(&mut self, node_ids: impl IntoIterator<Item = impl Borrow<usize>>) {
        for nid in node_ids {
            let nid = nid.borrow();
            // remove normal edges
            let biedge_keys: TinyVec<[EdgeKey;10]> = self.nodes
                .get_mut(nid).unwrap()
                .edges.iter().cloned()
                .collect();
            self.remove_biedges_from(&biedge_keys);
            // remove transitive edges
            let tedge_keys: TinyVec<[EdgeKey;10]> = self.nodes
                .get_mut(nid).unwrap()
                .transitives.iter().cloned()
                .collect();
            self.remove_transitives_from(&tedge_keys);
            self.nodes.remove(nid);
        }
    }

    fn remove_biedge(&mut self, edge_key:&EdgeKey) {
        let edge_key = biedge::canonical_edgekey(edge_key);
        self.edges.remove(&edge_key);
    }

    fn remove_biedges_from(&mut self, edge_keys: impl IntoIterator<Item = impl Borrow<EdgeKey>>) {
        for ekey in edge_keys {
            let ekey = ekey.borrow();
            if let Some(node) = self.nodes.get_mut(&ekey.id_from) {
                node.edges.retain(|key| key != ekey);
            }
            let ekey = biedge::flip_edgekey(ekey);
            if let Some(node) = self.nodes.get_mut(&ekey.id_from) {
                node.edges.retain(|key| key != &ekey);
            }
            self.remove_biedge(&ekey);
        }
    }

    fn remove_transitive(&mut self, edge_key:&EdgeKey) {
        let edge_key = biedge::canonical_edgekey(edge_key);
        self.transitives.remove(&edge_key);
    }

    fn remove_transitives_from(&mut self, edge_keys: impl IntoIterator<Item = impl Borrow<EdgeKey>>) {
        for ekey in edge_keys {
            let ekey = ekey.borrow();
            if let Some(node) = self.nodes.get_mut(&ekey.id_from) {
                node.transitives.retain(|key| key != ekey);
            }
            let ekey = biedge::flip_edgekey(ekey);
            if let Some(node) = self.nodes.get_mut(&ekey.id_from) {
                node.transitives.retain(|key| key != &ekey);
            }
            self.remove_transitive(&ekey);
        }
    }

    fn node_degree(&self, node_id:usize, dir:u8) -> usize {
        self.nodes[&node_id].edges.iter()
            .filter(|&key| key.strand_from == dir)
            .count()
    }

    fn edges_from(&self, node_id:usize, dir:u8) -> impl Iterator<Item = &EdgeKey> {
        self.nodes[&node_id].edges.iter()
            .filter(move |&key| key.strand_from == dir)
    }

    fn successors(&self, node_id:usize, dir:u8) -> impl Iterator<Item = (usize,u8)> {
        self.edges_from(node_id, dir)
            .map(|edge_key| (edge_key.id_to, edge_key.strand_to))
    }

    pub fn remove_weak_edges(&mut self, min_reads:usize) {

        let weak_edges = self.edges.values()
            .filter_map(|&edge_id| {
                let edge = self.get_biedge_idx(edge_id);
                if edge.nb_reads < min_reads {
                    // return Some(&edge.key)
                    // delete only edges that would not disconnect another node
                    let a_degree = self.node_degree(edge.key.id_from, edge.key.strand_from);
                    let b_degree = self.node_degree(edge.key.id_to, edge.key.strand_to);
                    if a_degree != 1 && b_degree != 1  {
                        return Some(&edge.key)
                    }
                }
                None
            }).cloned().collect_vec();
        self.remove_biedges_from(&weak_edges);

        let weak_bridges = self.edges.keys()
            .filter(|&edge_key| self.contains_transitive(edge_key))
            .filter(|&edge_key| {
                let edge_len = self.get_biedge(edge_key).unwrap().gaps.first().unwrap_or(&0).abs();
                let transitive_len = self.get_transitive(edge_key).unwrap().gaps.first().unwrap_or(&0).abs();
                let diff = 1.0 - (edge_len.min(transitive_len) as f64 / edge_len.max(transitive_len) as f64);
                spdlog::trace!("possible weak edge {edge_key}: edge_len={edge_len}, trans_len={transitive_len}, diff={diff:.3} weak={:?}", diff < 0.3);
                let a_degree = self.node_degree(edge_key.id_from, edge_key.strand_from);
                let b_degree = self.node_degree(edge_key.id_to, edge_key.strand_to);
                diff < 0.3 && a_degree != 1 && b_degree != 1
            })
            .cloned()
            .collect_vec();
        
        weak_bridges.iter().for_each(|edge_key| {
            let nb_reads = self.get_biedge(edge_key).unwrap().nb_reads;
            let transitive = self.get_transitive_mut(edge_key).unwrap();
            transitive.nb_reads += nb_reads;
            self.remove_biedge(edge_key);
        });
    }

    pub fn clear_transitive_edges(&mut self) {
        self.nodes.values_mut().for_each(|node| node.transitives.clear());
        self.transitives.clear();
    }

    /* GRAPH SIMPLIFICATION METHODS */

    pub fn add_bridges(&mut self, aware_alignments:&HashMap<usize,Vec<AwareAlignment>>) {
        self.clear_transitive_edges();
        // add bridges from aware alignments
        for (_query_idx, alignments) in aware_alignments.iter() {
            let mut first_indices: Vec<usize> = vec![];
            let mut last_indices: Vec<usize> = vec![];
            // find bifurcation nodes
            for (idx,(a,b)) in alignments.iter().tuple_windows().enumerate() {
                let EdgeKey { id_from, strand_from, id_to, strand_to } = EdgeKey::from_alignments(a,b);
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
                    let a = &alignments[first_idx];
                    let b = &alignments[last_idx];
                    if first_ctg_id != last_ctg_id && !a.is_ambiguous && !b.is_ambiguous {
                        let mut tedge_key = EdgeKey::from_alignments(a, b);
                        let was_canonical = tedge_key.canonicalize();
                        let tedge = self.get_transitive_or_create(&tedge_key);
                        tedge.nb_reads += 1;
                        tedge.min_shared_snvs = tedge.min_shared_snvs.max(a.nb_shared_snvs.min(b.nb_shared_snvs));
                        tedge.gaps.push((b.query_beg as i32) - (a.query_end as i32));
                        if tedge.seq_desc.is_empty() && b.query_beg > a.query_end {
                            let interval = SeqInterval{ tid: a.query_idx, beg: a.query_end, end: b.query_beg };
                            let strand = if was_canonical { b'+' } else { b'-' };
                            let aware_contig = AwareContig::new(SeqType::Read, interval, strand, 0);
                            tedge.seq_desc.push(aware_contig);
                        }
                    }
                }
            }
        }

        let nb_bridges = self.transitives.len();
        spdlog::debug!("{nb_bridges} read bridges added");
    }

    // TODO: the current implementation add all possible (non-ambiguous) transitive edges
    // The idea is to add only those that span nodes with outdegree > 1
    // so that I could possibly improve contiguity later
    pub fn add_bridges_2(&mut self, aware_alignments:&HashMap<usize,Vec<AwareAlignment>>) -> usize {
        self.clear_transitive_edges();

        for alignments in aware_alignments.values() {
            for (a_idx, a) in alignments.iter().enumerate().filter(|(_,a)| !a.is_ambiguous) {
                for b in alignments[a_idx+1..].iter().skip(1).filter(|b| !b.is_ambiguous) {
                    let mut tedge_key = EdgeKey::from_alignments(a, b);
                    let was_canonical = tedge_key.canonicalize();
                    let tedge = self.get_transitive_or_create(&tedge_key);
                    tedge.nb_reads += 1;
                    tedge.gaps.push((b.query_beg as i32) - (a.query_end as i32));
                    if tedge.seq_desc.is_empty() && b.query_beg > a.query_end {
                        let interval = SeqInterval{ tid: a.query_idx, beg: a.query_end, end: b.query_beg };
                        let strand = if was_canonical { b'+' } else { b'-' };
                        let aware_contig = AwareContig::new(SeqType::Read, interval, strand, 0);
                        tedge.seq_desc.push(aware_contig);
                    }
                }
            }
        }

        self.transitives.len()
    }

    fn find_junctions(&self) -> Vec<Junction> {
        let mut visited: HashSet<(usize,u8)> = HashSet::new();
        let mut junctions: Vec<Junction> = Vec::new();
        for node_id in self.nodes.keys().cloned() {
            for node_dir in [b'+',b'-'] {
                if visited.contains(&(node_id,node_dir)) {
                    continue
                }
                let in_edges = self.nodes[&node_id].edges.iter()
                    .filter(|key| key.strand_from == node_dir)
                    .cloned().collect_vec();
                if in_edges.len() <= 1 {
                    continue
                }
                // if !in_edges.iter().all(|edge_key| self.node_degree(edge_key.id_to, edge_key.strand_to) == 1) {
                //     continue
                // }

                // possibly found the "start" of junction
                visited.insert((node_id,node_dir));
                let mut junc = Junction::new();
                junc.in_edges.extend(in_edges);
                let mut in_dir = node_dir;
                let mut node = &self.nodes[&node_id];
                loop {
                    let out_dir = crate::seq::flip_strand(in_dir);
                    let out_edges = node.edges.iter()
                        .filter(|key| key.strand_from == out_dir)
                        .cloned().collect_vec();
                    if out_edges.is_empty() {
                        break;
                    }
                    let EdgeKey { id_from:src, strand_from:src_dir, id_to:succ, strand_to:succ_dir } = unsafe {
                        *out_edges.first().unwrap_unchecked()
                    };
                    if out_edges.len() > 1 { // found the "end" of the junction
                        visited.insert((src,src_dir));
                        junc.out_edges.extend(out_edges);

                        // if !junc.outputs().all(|(node_id,dir)| self.node_degree(node_id, dir) == 1) {
                        //     break;
                        // }

                        if HashSet::from_iter(junc.input_nodes()).intersection(&HashSet::from_iter(junc.output_nodes())).count() == 0 {
                            junctions.push(junc);
                        }

                        break;
                    }
                    node = &self.nodes[&succ];
                    if node.edges.iter().filter(|key| key.strand_from == succ_dir).count() > 1 {
                        break;
                    }
                    junc.mid_edges.push(out_edges[0]);
                    in_dir = succ_dir;
                };
            }
        }

        junctions
    }

    pub fn resolve_junctions(&mut self, min_reads:usize) {

        let mut tot_resolved = 0;
        for num_iter in 1_usize.. {

            spdlog::debug!("graph simplification iteration #{num_iter}");

            let mut junctions = self.find_junctions();
            junctions.sort_by_key(|j| j.inner_nodes().map(|x| self.nodes[&x].ctg.length()).sum::<usize>());

            let nb_resolved = junctions.iter().fold(0, |n,junction| {
                n + self.resolve_junction(junction, min_reads, num_iter) as usize
            });

            spdlog::debug!("junctions resolved: {nb_resolved}");

            if nb_resolved == 0 {
                break
            }

            tot_resolved += nb_resolved;
        }

        spdlog::debug!("total junctions resolved: {tot_resolved}");
    }

    fn resolve_junction(&mut self, junc:&Junction, min_reads:usize, num_iter:usize) -> bool {

        // resolve only junctions between phased sequences
        // if junc.inout_nodes().any(|node_id| !self.nodes[&node_id].ctg.is_phased()) {
        //     return false
        // }

        // first check if it is still a valid junction
        if junc.inputs().any(|(node_id,_)| !self.nodes.contains_key(&node_id)) || junc.outputs().any(|(node_id,_)| !self.nodes.contains_key(&node_id)) {
            return false
        }

        let mut bridges = vec![];
        for (node_id, node_dir) in junc.inputs() {
            for t_key in &self.nodes[&node_id].transitives {
                if (t_key.strand_from == node_dir) && junc.outputs().any(|out| out == (t_key.id_to,t_key.strand_to)) {
                    bridges.push(*t_key);
                }
            }
        }

        if num_iter == 1 {
            spdlog::trace!("{junc} with bridges:");
            for key in &bridges {
                let edge = self.get_transitive(key).unwrap();
                spdlog::trace!("  * {key} => nb_reads={}", edge.nb_reads);
            }
        }

        let ndeg: HashMap<usize, usize> = crate::utils::counter_from_iter(bridges.iter().flat_map(|key| [key.id_from,key.id_to]));

        let is_fully_covered = junc.inout_nodes().all(|n| ndeg.get(&n).is_some_and(|c| *c > 0));
        // let is_strictly_covered = is_fully_covered && bridges.iter().all(|key| {
        //     ndeg.get(&key.id_from).is_some_and(|c| *c == 1) || ndeg.get(&key.id_to).is_some_and(|c| *c == 1)
        // });

        if !is_fully_covered {
            // if !bridges.is_empty() {
            //     spdlog::debug!("Uncovered junction {junc}");
            //     for key in &bridges {
            //         let bridge = self.get_transitive(key).unwrap();
            //         spdlog::debug!("  * {key} (reads: {}, SNVs:{})", bridge.nb_reads, bridge.min_shared_snvs);
            //     }
            // }
            self.remove_transitives_from(&bridges);
            return false
        }

        // if is_fully_covered && !is_strictly_covered {
        //     if !bridges.is_empty() {
        //         spdlog::debug!("Junction not strictly covered: {junc}");
        //         for key in &bridges {
        //             let bridge = self.get_transitive(key).unwrap();
        //             spdlog::debug!("  * {key} (reads: {}, SNVs:{})", bridge.nb_reads, bridge.min_shared_snvs);
        //         }
        //     }
        // }

        // remove low-weight edges while keeping the junction fully covered
        bridges.sort_by_key(|key| -(self.get_transitive(key).unwrap().nb_reads as i32));
        loop {
            if self.get_transitive(bridges.last().unwrap()).unwrap().nb_reads >= min_reads {
                break
            }
            let ndeg: HashMap<usize, usize> = crate::utils::counter_from_iter(bridges[..bridges.len()-1].iter().flat_map(|key| [key.id_from,key.id_to]));
            let is_fully_covered = junc.inout_nodes().all(|n| ndeg.get(&n).is_some_and(|c| *c > 0));
            if is_fully_covered {
                bridges.pop();
            } else {
                break
            }
        }

        // if is_fully_covered && !is_strictly_covered {
        //     if !bridges.is_empty() {
        //         spdlog::debug!("Solving with:");
        //         for key in &bridges {
        //             let bridge = self.get_transitive(key).unwrap();
        //             spdlog::debug!("  * {key} (reads: {}, SNVs:{})", bridge.nb_reads, bridge.min_shared_snvs);
        //         }
        //     }
        // }

        // if !is_strictly_covered {
        //     // try remove low-weight edges and see if it becomes strictly covered.
        //     bridges.retain(|key| self.get_transitive(key).is_some_and(|e| e.observations >= min_reads));
        //     let ndeg: HashMap<usize, usize> = crate::utils::counter_from_iter(bridges.iter().flat_map(|key| [key.id_from,key.id_to]));
        //     if junc.inout_nodes().any(|n| ndeg.get(&n).is_none_or(|c| *c == 0)) {
        //         self.remove_transitives_from(&bridges);
        //         return false
        //     }
        // }

        let mut bridge_paths = vec![];
        for EdgeKey { id_from, strand_from, id_to, strand_to } in &bridges {
            let in_edge = biedge::flip_edgekey(junc.in_edges.iter().find(|key| (&key.id_to,&key.strand_to) == (id_from,strand_from)).unwrap());
            let mut path = Vec::with_capacity(junc.mid_edges.len()+2);
            path.push(in_edge);
            path.extend(junc.mid_edges.iter());
            let out_edge = *junc.out_edges.iter().find(|key| (&key.id_to,&key.strand_to) == (id_to,strand_to)).unwrap();
            path.push(out_edge);
            bridge_paths.push(path);
        }

        for path in &bridge_paths {
            assert!(!path.is_empty());
            let first = unsafe { path.first().unwrap_unchecked() };
            let last = unsafe { path.last().unwrap_unchecked() };
            let edge_key = EdgeKey::new(first.id_from, first.strand_from, last.id_to, last.strand_to);
            if self.contains_edge(&edge_key) {
                // eprintln!("edge {edge_key} already exists for {junc}");
                return false
            }
            if path.iter().any(|edge_key| !self.contains_edge(edge_key)) { // path does not exist
                return false
            }
        }

        assert!(!bridge_paths.is_empty());
        self.resolve_bridge_paths(&bridge_paths);
        self.remove_transitives_from(&bridges);
        true
    }

    fn resolve_bridge_paths(&mut self, bridge_paths:&[Vec<EdgeKey>]) {
        
        let deleted_nodes: HashSet<usize> = bridge_paths.iter()
            .flat_map(|path| path.iter().skip(1))
            .map(|key| key.id_from)
            .collect();
        let deleted_edges: HashSet<EdgeKey> = bridge_paths.iter().flatten().cloned().collect();

        for path in bridge_paths {
            assert!(path.len() >= 2);
            let first = unsafe { path.first().unwrap_unchecked() };
            let last = unsafe { path.last().unwrap_unchecked() };
            let tkey = EdgeKey::new(first.id_from, first.strand_from, last.id_to, last.strand_to);
            let tedge = self.get_transitive(&tkey).unwrap();
            let tedge_nb_reads = tedge.nb_reads;
            let tedge_min_shared_snvs = tedge.min_shared_snvs;
            let tedge_gaps = tedge.gaps.clone();
            let tedge_seq_desc = tedge.seq_desc.clone();

            // create new edge
            let edge = self.get_biedge_or_create(&tkey);
            edge.nb_reads = tedge_nb_reads;
            edge.min_shared_snvs = tedge_min_shared_snvs;
            edge.gaps = tedge_gaps; // junc_seq
            edge.seq_desc = tedge_seq_desc;
        }
        // delete old nodes/edges
        self.remove_biedges_from(&deleted_edges);
        self.remove_nodes_from(&deleted_nodes);
    }

    fn expand_edges(&mut self) {
        
        let nonempty_edges = self.edges.iter()
            .filter_map(|(key,id)| if !self.edge_data[*id].seq_desc.is_empty() { Some(key) } else { None })
            .cloned().collect_vec();

        for edge_key in nonempty_edges.iter() {
            let edge = self.get_biedge(edge_key).unwrap();
            assert!(edge_key.is_canonical());
            let node_id = self.add_node(edge.seq_desc[0]);
            let new_edge_key = EdgeKey::new(edge_key.id_from, edge_key.strand_from, node_id, b'+');
            self.get_biedge_or_create(&new_edge_key);
            let new_edge_key = EdgeKey::new( node_id, b'-', edge_key.id_to, edge_key.strand_to);
            self.get_biedge_or_create(&new_edge_key);
        }

        self.remove_biedges_from(nonempty_edges);
    }


    fn extract_aware_paths(&self) -> (Vec<AwarePath>, HashMap<usize,usize>) {
        let mut visited: HashSet<usize> = HashSet::new();
        let mut aware_paths: Vec<AwarePath> = Vec::new();
        let mut node_to_path: HashMap<usize, usize> = HashMap::new();
        for node_id in self.nodes.keys() {
            if visited.insert(*node_id) {
                let edges = self.unitig_from(*node_id, &mut visited);
                let nodes = if edges.is_empty() {
                    vec![*node_id]
                } else {
                    let first = unsafe { edges.first().unwrap_unchecked().id_from };
                    std::iter::once(first).chain(edges.iter().map(|key| key.id_to)).collect_vec()
                };
                let id = aware_paths.len();
                nodes.iter().for_each(|n| { node_to_path.insert(*n, id); });
                aware_paths.push(AwarePath{ nodes, edges });
            }
        }
        (aware_paths, node_to_path)
    }

    fn unitig_from(&self, start_id:usize, visited: &mut HashSet<usize>) -> Vec<EdgeKey> {

        let mut linear_path_from = |mut node_id:usize, mut node_dir:u8| -> Vec<EdgeKey> {
            let mut path_edges = vec![];
            loop {
                let node = &self.nodes[&node_id];
                let mut outedges = node.edges.iter().cloned().filter(|key| key.strand_from == node_dir).take(2).collect::<ArrayVec<[EdgeKey;2]>>();
                if outedges.len() != 1 {
                    break
                }
                let outkey = unsafe { outedges.pop().unwrap_unchecked() };
                if visited.contains(&outkey.id_to) || self.nodes[&outkey.id_to].edges.iter().filter(|key| key.strand_from == outkey.strand_to).count() > 1 {
                    break
                }
                (node_id, node_dir) = (outkey.id_to, crate::seq::flip_strand(outkey.strand_to));
                path_edges.push(outkey);
                visited.insert(node_id);
            }
            path_edges
        };
        
        let mut path_edges = linear_path_from(start_id, b'+');
        path_edges.reverse();
        path_edges.iter_mut().for_each(|key| key.flip());
        path_edges.append(&mut linear_path_from(start_id, b'-'));
        path_edges
    }


    pub fn build_haplotigs(&self, aware_paths: &[AwarePath], ref_db: &SeqDatabase, read_db: &SeqDatabase, fragment_dir: &Path, work_dir: &Path, opts: &Options) -> Result<Vec<BitSeq>> {

        let mut unpolished_bp: usize = 0;
        let mut haplotigs = Vec::with_capacity(aware_paths.len());
        for (id, aware_path) in aware_paths.iter().enumerate() {

            assert!(!aware_path.nodes.is_empty());

            // create backbone sequence
            
            let backbone_sequence = {

                let mut backbone_sequence = Vec::new();
                let mut is_phased = false;

                let node_id = unsafe { aware_path.nodes.first().unwrap_unchecked() };
                let node_dir = if aware_path.edges.is_empty() { b'+' } else { crate::seq::flip_strand(aware_path.edges[0].strand_from) };
                let node_iter = std::iter::once((*node_id, node_dir))
                    .chain(aware_path.edges.iter().map(|key| (key.id_to, key.strand_to)));

                // TODO: handle negative gaps between sequences?
                for (node_id, node_dir) in node_iter {
                    let aware_contig = &self.nodes[&node_id].ctg;
                    let SeqInterval { tid, beg, end } = aware_contig.interval();
                    let mut sequence = match aware_contig.contig_type() {
                        SeqType::Haplotype(_) | SeqType::Unphased => { ref_db.sequences[tid][beg..end].to_vec() },
                        SeqType::Read => { read_db.sequences[tid][beg..end].to_vec() },
                        // SeqType::Haplotype(_) | SeqType::Unphased => { ref_db.sequences[tid].subseq(beg, end) },
                        // SeqType::Read => { read_db.sequences[tid].subseq(beg, end) },
                    };
                    if node_dir != aware_contig.strand() { crate::seq::revcomp_inplace(&mut sequence); }
                    is_phased = is_phased || aware_contig.is_phased();
                    backbone_sequence.append(&mut sequence);
                }

                // if path does not contain phased regions, do not polish
                if aware_path.nodes.iter().all(|id| ! matches!(self.nodes[id].ctg.contig_type(), SeqType::Haplotype(_)) ) {
                    haplotigs.push(BitSeq::from_utf8(&backbone_sequence));
                    continue
                }

                backbone_sequence
            };

            // create specific directory for polishing current sequence

            let tmp_dir = work_dir.join(format!("{}_tmp", id));
            std::fs::create_dir_all(&tmp_dir)
                .with_context(|| format!("Cannot create output directory: \"{}\"", tmp_dir.display()))?;

            // output backbone sequence

            let target_path = {
                let out_path = tmp_dir.join(format!("{}.unpolished.fa.gz", id));
                let mut writer = crate::utils::get_file_writer(&out_path);
                writer.write_all(format!(">ctg{id}\n").as_bytes())?;
                writer.write_all(&backbone_sequence)?;
                writer.write_all(b"\n")?;
                out_path
            };
            
            // output haplotype-specific reads associated to the target sequence

            let read_path = {
                let out_path = tmp_dir.join(format!("{}.reads.fa.gz", id));
                let mut read_ids = HashSet::new();
                let mut unitig_read_writer = crate::utils::get_file_writer(&out_path);
                let haplotype_files = aware_path.nodes.iter()
                    .filter_map(|nid| self.nodes[nid].ctg.haplotype_id())
                    .map(|h| fragment_dir.join(format!("{}_{}-{}_h{}.fa.gz", ref_db.names[h.tid], h.beg, h.end, h.hid)));
                for htfile in haplotype_files {
                    if let Ok(mut reader) = needletail::parse_fastx_file(&htfile) {
                        while let Some(record) = reader.next() {
                            let record = record.unwrap();
                            let read_name = std::str::from_utf8(record.id().split(|b| b.is_ascii_whitespace()).next().unwrap()).unwrap();
                            if !read_ids.contains(read_name) {
                                read_ids.insert(read_name.to_string());
                                record.write(&mut unitig_read_writer, None)?;
                            }
                        }
                    } else {
                        spdlog::warn!("Empty haplotype read file: {}", htfile.display());
                    }
                }
                out_path
            };

            // run racon to polish the backbone sequence with strain-specific reads

            let polished_path = tmp_dir.join(format!("{}.polished.fa", id));
            crate::polish::racon_polish(&target_path, &read_path, &polished_path, PolishMode::Oblivious, &tmp_dir, opts)?;
            
            let polished_sequence = {
                // let mut reader = needletail::parse_fastx_file(&polished_path)
                //     .with_context(||format!("Cannot read polished sequence from: {}", polished_path.display()))?;
                match needletail::parse_fastx_file(&polished_path) {
                    Ok(mut reader) => {
                        if let Some(Ok(record)) = reader.next() {
                            BitSeq::from_utf8(&record.seq())
                        } else {
                            spdlog::trace!("unitig {id} was not polished ({} bp)", backbone_sequence.len());
                            unpolished_bp += backbone_sequence.len();
                            BitSeq::from_utf8(&backbone_sequence)
                        }
                    },
                    Err(_) => { BitSeq::from_utf8(&backbone_sequence) }
                }
            };

            haplotigs.push(polished_sequence);

            if !opts.keep_temp {
                std::fs::remove_dir_all(&tmp_dir)?;
            }
        }

        spdlog::debug!("unpolished bps: {unpolished_bp}");

        Ok(haplotigs)
    }

    fn write_info(&self, ref_db: &SeqDatabase, read_db: &SeqDatabase, aware_paths:&[AwarePath], info_path:&Path) -> Result<()> {
        
        let mut info_writer = crate::utils::get_file_writer(info_path);
        
        for (utg_id, ap) in aware_paths.iter().enumerate() {
            let node_id = unsafe { ap.nodes.first().unwrap_unchecked() };
            let node_dir = if ap.edges.is_empty() { b'+' } else { crate::seq::flip_strand(ap.edges[0].strand_from) };
            let node_iter = std::iter::once((*node_id, node_dir))
                .chain(ap.edges.iter().map(|key| (key.id_to, key.strand_to)));

            let mut utg_pos = 0;
            for (node_id, node_dir) in node_iter {

                let aware_contig = &self.nodes[&node_id].ctg;
                let SeqInterval { tid, beg, end } = aware_contig.interval();
                let is_reverse = node_dir != aware_contig.strand();
                let strand = b"+-"[is_reverse as usize] as char;
                let ctg_len = end - beg;
                let utg_end = utg_pos + ctg_len;

                match aware_contig.contig_type() {
                    SeqType::Unphased => {
                        let ref_name = ref_db.names[tid].as_str();
                        let line = format!("{utg_id}\t{utg_pos}\t{utg_end}\t{ctg_len}\tREF\t{ref_name}\t{beg}\t{end}\t{strand}\n");
                        info_writer.write_all(line.as_bytes())?;
                    },
                    SeqType::Haplotype(hid) => {
                        let ref_name = ref_db.names[tid].as_str();
                        let line = format!("{utg_id}\t{utg_pos}\t{utg_end}\t{ctg_len}\tHAP\t{ref_name}\t{beg}\t{end}\th{hid}\t{strand}\n");
                        info_writer.write_all(line.as_bytes())?;
                    },
                    SeqType::Read => {
                        let read_name = read_db.names[tid].as_str();
                        let line = format!("{utg_id}\t{utg_pos}\t{utg_end}\t{ctg_len}\tREAD\t{read_name}\t{beg}\t{end}\t{strand}\n");
                        info_writer.write_all(line.as_bytes())?;
                    },
                };

                utg_pos = utg_end;
            }
        }

        Ok(())
    }

    // TODO: handle self loops
    // TODO: handle negative gaps
    pub fn build_assembly_graph(&mut self, ref_db: &SeqDatabase, read_db: &SeqDatabase, fragment_dir:&Path, work_dir:&Path, opts: &Options) -> Result<AsmGraph<usize>> {

        std::fs::create_dir_all(work_dir)
            .with_context(|| format!("Cannot create output directory: \"{}\"", work_dir.display()))?;

        self.expand_edges();

        let (aware_paths, node_index) = self.extract_aware_paths();

        let info_path = work_dir.join("paths.info.txt");
        self.write_info(ref_db, read_db, &aware_paths, &info_path)?;
        
        let haplotigs = self.build_haplotigs(&aware_paths, ref_db, read_db, fragment_dir, work_dir, opts)?;
        assert!(haplotigs.len() == aware_paths.len());

        // build assembly graph 

        let mut asm_graph = super::asmgraph::AsmGraph::new();
        for (id, sequence) in haplotigs.into_iter().enumerate() {
            asm_graph.add_node(id, sequence);
        }

        // TODO: handle self loops
        for edge_key in self.edges.keys() {

            // if !node_index.contains_key(&edge_key.id_from) || !node_index.contains_key(&edge_key.id_to) {
            //     continue
            // }
            
            let id_from = node_index[&edge_key.id_from];
            let id_to = node_index[&edge_key.id_to];

            if id_from != id_to || edge_key.id_from == edge_key.id_to {
                let path_from = &aware_paths[id_from];
                let strand_from = if path_from.nodes.len() == 1 {
                    edge_key.strand_from
                } else {
                    b"-+"[(edge_key.id_from == path_from.nodes[0]) as usize]
                };
                let path_to = &aware_paths[id_to];
                let strand_to = if path_to.nodes.len() == 1 {
                    edge_key.strand_to
                } else {
                    b"-+"[(edge_key.id_to == path_to.nodes[0]) as usize]
                };
                asm_graph.add_link(Link::new(id_from, strand_from, id_to, strand_to));
            }
        }

        Ok(asm_graph)
    }


    /* OUTPUT METHODS */

    pub fn write_gfa(&self, gfa_path:PathBuf, ref_db:&SeqDatabase) -> std::io::Result<()> {
        let mut gfa = crate::utils::get_file_writer(gfa_path.as_path());
        gfa.write_all(b"H\tVN:Z:1.0\n")?;
        for (node_id, node) in self.nodes.iter() {
            let node_ctg = node.ctg;
            let node_name = match node_ctg.contig_type() {
                SeqType::Haplotype(hid) => format!("{}_{}-{}_h{}_id{}", ref_db.names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), hid, node_id),
                SeqType::Unphased => format!("{}_{}-{}_id{}", ref_db.names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), node_id),
                SeqType::Read => format!("read{}_{}-{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), node_id),
            };
            let node_line = format!("S\t{}\t*\tLN:i:{}\tdp:f:{:.1}\n", node_name, node_ctg.length(), node_ctg.depth());
            gfa.write_all(node_line.as_bytes())?;
        }

        for edge_id in self.edges.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let node_ctg = self.nodes[&id_from].ctg;
            let name_from = match node_ctg.contig_type() {
                SeqType::Haplotype(hid) => format!("{}_{}-{}_h{}_id{}", ref_db.names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), hid, id_from),
                SeqType::Unphased => format!("{}_{}-{}_id{}", ref_db.names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), id_from),
                SeqType::Read => format!("read{}_{}-{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), id_from),
            };
            let node_ctg = self.nodes[&id_to].ctg;
            let name_to = match node_ctg.contig_type() {
                SeqType::Haplotype(hid) => format!("{}_{}-{}_h{}_id{}", ref_db.names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), hid, id_to),
                SeqType::Unphased => format!("{}_{}-{}_id{}", ref_db.names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), id_to),
                SeqType::Read => format!("read{}_{}-{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), id_to),
            };
            let edge_line = format!("L\t{}\t{}\t{}\t{}\t0M\tRC:i:{}\n", name_from, crate::seq::flip_strand(strand_from) as char, name_to, strand_to as char, edge.nb_reads );
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
            let edge_label = format!("{}/{}bp", edge.nb_reads, edge_gap);
            let edge_line = format!("\t{id_from} -> {id_to}\t[arrowtail={arrow_tail}, arrowhead={arrow_head}, dir=both, label=\"{edge_label}\"];\n");
            dot.write_all(edge_line.as_bytes())?;
        }
        for edge_id in self.transitives.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let arrow_tail = if strand_from == b'+' { "normal" } else { "inv" };
            let arrow_head = if strand_to == b'+' { "normal" } else { "inv" };
            let edge_gap = if edge.gaps.is_empty() { 0 } else { *edge.gaps.iter().sorted_unstable().nth(edge.gaps.len()/2).unwrap() };
            let edge_label = format!("{}/{}bp", edge.nb_reads, edge_gap);
            let edge_line = format!("\t{id_from} -> {id_to}\t[arrowtail={arrow_tail}, arrowhead={arrow_head}, dir=both, color=\"red\", label=\"{edge_label}\"];\n");
            dot.write_all(edge_line.as_bytes())?;
        }
        dot.write_all(b"}\n")
    }

}

