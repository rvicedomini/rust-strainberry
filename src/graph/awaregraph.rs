use std::borrow::Borrow;
use std::path::PathBuf;
use itertools::Itertools;
use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use tinyvec::{TinyVec, ArrayVec};

use crate::awarecontig::{AwareAlignment, AwareContig, SeqType};
use crate::graph::asmgraph::AsmGraph;
use crate::seq::SeqInterval;
use crate::utils::flip_strand;

use super::biedge::{self, BiEdge, EdgeKey, Node};
use super::junction::Junction;

#[derive(Debug)]
pub struct Unitig {
    pub id: usize,
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
                let canon_key = biedge::canonical_edgekey(&edge_key);
                let edge = aware_graph.get_biedge_or_create(&canon_key);
                let interval = SeqInterval{ tid: a.tid(), beg: a.end(), end: b.beg() };
                // strand of edge sequence is defined according to the "direction" of the canonical edge
                let aware_contig = AwareContig::new(SeqType::Unphased, interval, b'-', 0);
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
            for (a, b) in alignments.iter().tuple_windows() {
                self.add_edge_from_consecutive_alignments(a,b);
            }
        }
    }

    fn add_edge_from_consecutive_alignments(&mut self, a:&AwareAlignment, b:&AwareAlignment) {
        let mut edge_key: EdgeKey = EdgeKey::from_alignments(a, b);
        let was_canonical = edge_key.canonicalize();
        let edge = self.get_biedge_or_create(&edge_key);
        edge.observations += 1;
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

    pub fn clear_transitive_edges(&mut self) {
        self.nodes.values_mut().for_each(|node| node.transitives.clear());
        self.transitives.clear();
    }

    /* GRAPH SIMPLIFICATION METHODS */

    pub fn add_bridges(&mut self, aware_alignments:&HashMap<usize,Vec<AwareAlignment>>) -> usize {
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
                    if first_ctg_id != last_ctg_id {
                        let a = &alignments[first_idx];
                        let b = &alignments[last_idx];
                        let mut tedge_key = EdgeKey::from_alignments(a, b);
                        let was_canonical = tedge_key.canonicalize();
                        let tedge = self.get_transitive_or_create(&tedge_key);
                        tedge.observations += 1;
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
                // possibly found the "start" of junction
                visited.insert((node_id,node_dir));
                let mut junc = Junction::new();
                junc.in_edges.extend(in_edges);
                let mut in_dir = node_dir;
                let mut node = &self.nodes[&node_id];
                loop {
                    let out_dir = flip_strand(in_dir);
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

    pub fn resolve_read_bridges(&mut self, min_reads:usize) -> usize {
        let mut resolved = 0;
        let mut junctions = self.find_junctions();
        junctions.sort_unstable_by_key(|j| j.size());
        for junc in &junctions {
            // println!("Resolving junction: {junc}");
            resolved += self.resolve_read_junction(junc, min_reads) as usize;
        }
        // println!("Junction resolution: found={} resolved={resolved}", junctions.len());
        resolved
    }

    fn resolve_read_junction(&mut self, junc:&Junction, min_reads:usize) -> bool {
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

        let ndeg: HashMap<usize, usize> = crate::utils::counter_from_iter(bridges.iter().flat_map(|key| [key.id_from,key.id_to]));

        let is_fully_covered = junc.inout_nodes().all(|n| ndeg.get(&n).is_some_and(|c| *c > 0));
        let is_strictly_covered = is_fully_covered && bridges.iter().all(|key| {
            ndeg.get(&key.id_from).is_some_and(|c| *c == 1) || ndeg.get(&key.id_to).is_some_and(|c| *c == 1)
        });

        if !is_fully_covered {
            // println!("  not fully covered");
            self.remove_transitives_from(&bridges);
            return false
        }

        if !is_strictly_covered {
            // println!("  not strictly covered");
            bridges.retain(|key| self.get_transitive(key).is_some_and(|e| e.observations >= min_reads));
            let ndeg: HashMap<usize, usize> = crate::utils::counter_from_iter(bridges.iter().flat_map(|key| [key.id_from,key.id_to]));
            if junc.inout_nodes().any(|n| ndeg.get(&n).is_none_or(|c| *c == 0)) {
                self.remove_transitives_from(&bridges);
                return false
            }
        }

        // println!("  strictly covered with bridges:");
        // for edge in &bridges {
        //     println!("    * {edge}");
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

        // for key in bridges:
        //     tedge = self._get_tedge(key)
        //     if (key[0],key[1]) not in tedge.gapseq and junc.size() > 1:
        //         # logger.debug(f'{junc} will not be resolved')
        //         return False

        self.resolve_read_bridged_paths(&bridge_paths,junc);
        self.remove_transitives_from(&bridges);

        true
    }

    fn resolve_read_bridged_paths(&mut self, bridge_paths:&[Vec<EdgeKey>], _junc:&Junction) {

        // println!("  resolving with:");
        // for path in bridge_paths {
        //     println!("    * {path:?}");
        // }
        
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
            let tedge_observations = tedge.observations;
            let tedge_gaps = tedge.gaps.clone();
            let tedge_seq_desc = tedge.seq_desc.clone();

            // create new edge
            let edge = self.get_biedge_or_create(&tkey);
            edge.observations = tedge_observations;
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


    fn unitig_partition(&self) -> (Vec<Unitig>, HashMap<usize,usize>) {
        let mut visited: HashSet<usize> = HashSet::new();
        let mut unitigs: Vec<Unitig> = Vec::new();
        let mut node_to_unitig: HashMap<usize, usize> = HashMap::new();
        for node_id in self.nodes.keys() {
            if visited.insert(*node_id) {
                let edges = self.unitig_from(*node_id, &mut visited);
                let nodes = if edges.is_empty() {
                    vec![*node_id]
                } else {
                    let first = unsafe { edges.first().unwrap_unchecked().id_from };
                    std::iter::once(first).chain(edges.iter().map(|key| key.id_to)).collect_vec()
                };
                let unitig_id = unitigs.len();
                nodes.iter().for_each(|n| { node_to_unitig.insert(*n, unitig_id); });
                unitigs.push(Unitig{ id:unitig_id, nodes, edges });
            }
        }
        (unitigs, node_to_unitig)
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
                (node_id, node_dir) = (outkey.id_to, flip_strand(outkey.strand_to));
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

    // TODO: handle self loops
    // TODO: handle negative gaps
    pub fn compact_graph(&mut self, target_sequences: &[Vec<u8>], read_sequences: &HashMap<usize,Vec<u8>>) -> super::asmgraph::AsmGraph<Unitig> {

        self.expand_edges();

        let (unitigs, node_to_unitig) = self.unitig_partition();
        let mut unitig_graph: AsmGraph<Unitig> = super::asmgraph::AsmGraph::new();

        for unitig in unitigs {
            // TODO: build unitig sequence from nodes/edges of the path
            if unitig.nodes.len() == 1 {
                let node_id= unsafe { unitig.nodes.first().unwrap_unchecked() };
                let aware_contig = &self.nodes[node_id].ctg;
                let SeqInterval { tid, beg, end } = aware_contig.interval();
                let mut sequence = match aware_contig.contig_type() {
                    SeqType::Haplotype(_) | SeqType::Unphased => { target_sequences[tid][beg..end].to_vec() },
                    SeqType::Read => { read_sequences[&tid][beg..end].to_vec() },
                };
                if aware_contig.strand() == b'-' { crate::seq::read::revcomp_inplace(&mut sequence); }
                unitig_graph.add_node(unitig.id, sequence, Some(unitig));
                continue
            }
            // TODO: handle negative gaps
            assert!(unitig.nodes.len() > 1);
            let mut unitig_sequence = vec![];
            let node_iter = std::iter::once(
                    (unitig.edges[0].id_from, flip_strand(unitig.edges[0].strand_from))
                ).chain(unitig.edges.iter().map(|key| (key.id_to, key.strand_to)));
            for (node_id, node_dir) in node_iter {
                let aware_contig = &self.nodes[&node_id].ctg;
                let SeqInterval { tid, beg, end } = aware_contig.interval();
                let mut sequence = match aware_contig.contig_type() {
                    SeqType::Haplotype(_) | SeqType::Unphased => { target_sequences[tid][beg..end].to_vec() },
                    SeqType::Read => { read_sequences[&tid][beg..end].to_vec() },
                };
                if node_dir != aware_contig.strand() { crate::seq::read::revcomp_inplace(&mut sequence); }
                unitig_sequence.append(&mut sequence);
            }
            unitig_graph.add_node(unitig.id, unitig_sequence, Some(unitig));
        }

        for edge_key in self.edges.keys() {
            let node_from = unitig_graph.get_node(node_to_unitig[&edge_key.id_from]).unwrap();
            let node_to = unitig_graph.get_node(node_to_unitig[&edge_key.id_to]).unwrap();
            let utg_from = node_from.data.as_ref().unwrap();
            let utg_to = node_to.data.as_ref().unwrap();
            if utg_from.id != utg_to.id || edge_key.id_from == edge_key.id_to {
                let dir_from = if utg_from.nodes.len() == 1 {
                    edge_key.strand_from
                } else {
                    b"-+"[(edge_key.id_from == utg_from.nodes[0]) as usize]
                };
                let dir_to = if utg_to.nodes.len() == 1 {
                    edge_key.strand_to
                } else {
                    b"-+"[(edge_key.id_to == utg_to.nodes[0]) as usize]
                };
                unitig_graph.add_edge(EdgeKey::new(utg_from.id, dir_from, utg_to.id, dir_to));
            }
        }

        unitig_graph
    }


    /* OUTPUT METHODS */

    pub fn write_gfa(&self, gfa_path:PathBuf, target_names:&[String]) -> std::io::Result<()> {
        let mut gfa = crate::utils::get_file_writer(gfa_path.as_path());
        gfa.write_all(b"H\tVN:Z:1.0\n")?;
        for (node_id, node) in self.nodes.iter() {
            let node_ctg = node.ctg;
            let node_name = match node_ctg.contig_type() {
                SeqType::Haplotype(hid) => format!("{}_{}-{}_h{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), hid, node_id),
                SeqType::Unphased => format!("{}_{}-{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), node_id),
                SeqType::Read => format!("read{}_{}-{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), node_id),
            };
            let node_line = format!("S\t{}\t*\tLN:i:{}\tdp:i:{}\n", node_name, node_ctg.length(), node_ctg.depth() as usize);
            gfa.write_all(node_line.as_bytes())?;
        }

        for edge_id in self.edges.values() {
            let edge = self.get_biedge_idx(*edge_id);
            let EdgeKey { id_from, strand_from, id_to, strand_to } = edge.key;
            let node_ctg = self.nodes[&id_from].ctg;
            let name_from = match node_ctg.contig_type() {
                SeqType::Haplotype(hid) => format!("{}_{}-{}_h{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), hid, id_from),
                SeqType::Unphased => format!("{}_{}-{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), id_from),
                SeqType::Read => format!("read{}_{}-{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), id_from),
            };
            let node_ctg = self.nodes[&id_to].ctg;
            let name_to = match node_ctg.contig_type() {
                SeqType::Haplotype(hid) => format!("{}_{}-{}_h{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), hid, id_to),
                SeqType::Unphased => format!("{}_{}-{}_id{}", target_names[node_ctg.tid()], node_ctg.beg(), node_ctg.end(), id_to),
                SeqType::Read => format!("read{}_{}-{}_id{}", node_ctg.tid(), node_ctg.beg(), node_ctg.end(), id_to),
            };
            let edge_line = format!("L\t{}\t{}\t{}\t{}\t0M\tRC:i:{}\n", name_from, crate::utils::flip_strand(strand_from) as char, name_to, strand_to as char, edge.observations );
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

