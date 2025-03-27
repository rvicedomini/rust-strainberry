use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use itertools::Itertools;
use ahash::AHashMap as HashMap;
use ahash::AHashSet as HashSet;
use tinyvec::{tiny_vec, TinyVec};

use crate::utils::BamRecordId;
use super::haplotype::{HaplotypeId, Haplotype};

#[derive(Debug, Default, Clone)]
struct Edge(HaplotypeId,usize);

type AdjList = HashMap<HaplotypeId,TinyVec<[Edge;10]>>;

pub struct HaploGraph {
    haplotypes: HashMap<HaplotypeId,Haplotype>,
    succ: AdjList,
    pred: AdjList,
}

impl HaploGraph {

    //TODO: do not allow "transitive" edges (due to reads not assigned to haplotypes when ambiguous)
    pub fn new(haplotypes: HashMap<HaplotypeId,Haplotype>, sread_haplotypes: &HashMap<BamRecordId, Vec<HaplotypeId>>) -> HaploGraph {
        let mut succ: AdjList = AdjList::default();
        let mut pred: AdjList = AdjList::default();
        
        for ht_uid in haplotypes.keys() {
            succ.insert(*ht_uid, tiny_vec!());
            pred.insert(*ht_uid, tiny_vec!());
        }

        for ht_list in sread_haplotypes.values() {
            for (a,b) in ht_list.iter().sorted_unstable_by_key(|&hid| haplotypes[hid].beg()).tuple_windows() {
                let a_succ = succ.get_mut(a).unwrap();
                if let Some(idx) = a_succ.iter().position(|Edge(uid,_)| uid == b) {
                    a_succ[idx].1 += 1;
                } else {
                    a_succ.push(Edge(*b,1));
                }
                let b_pred = pred.get_mut(b).unwrap();
                if let Some(idx) = b_pred.iter().position(|Edge(uid,_)| uid == a) {
                    b_pred[idx].1 += 1;
                } else {
                    b_pred.push(Edge(*a,1));
                }
            }
        }

        HaploGraph{ haplotypes, succ, pred }
    }

    fn out_degree(&self, hid:&HaplotypeId) -> usize { self.succ[hid].len() }
    fn in_degree(&self, hid:&HaplotypeId) -> usize { self.pred[hid].len() }

    pub fn remove_weak_edges(&mut self, min_reads:usize, min_frac:f64) {

        let mut deleted_edges = HashSet::new();
        for node in self.succ.keys() {
            if self.out_degree(node) > 0 {
                let out_reads: usize = self.succ[node].iter().map(|&Edge(_,n)| n).sum();
                self.succ[node].iter()
                    .filter(|Edge(_,edge_reads)| *edge_reads < min_reads || (*edge_reads as f64)/(out_reads as f64) < min_frac)
                    .for_each(|Edge(succ_id,_)| {
                        deleted_edges.insert((*node,*succ_id));
                    });
            }
            if self.in_degree(node) > 0 {
                let in_reads: usize = self.pred[node].iter().map(|&Edge(_,n)| n).sum();
                self.pred[node].iter()
                    .filter(|Edge(_,edge_reads)| *edge_reads < min_reads || (*edge_reads as f64)/(in_reads as f64) < min_frac)
                    .for_each(|Edge(pred_id,_)| {
                        deleted_edges.insert((*pred_id,*node));
                    });
            }
        }
        deleted_edges.drain().for_each(|(a,b)| self.delete_edge(a,b));

        let mut phased_map: HashMap<(usize,usize),TinyVec<[HaplotypeId;10]>> = HashMap::new();
        for hid in self.succ.keys() {
            phased_map.entry((hid.beg, hid.end))
                .or_default()
                .push(*hid);
        }

        let sorted_regions = phased_map.keys().cloned().sorted();
        for (prev,succ) in sorted_regions.tuple_windows() {
            let nb_edges: usize = phased_map[&prev].iter().map(|hid| self.out_degree(hid)).sum();
            let prev_nodes: usize = phased_map[&prev].len();
            let succ_nodes: usize = phased_map[&succ].len();
            if nb_edges != succ_nodes || prev_nodes != succ_nodes {
                deleted_edges.extend(phased_map[&prev].iter()
                    .flat_map(|a| {
                        self.succ.get(a).unwrap().iter()
                            .map(move |Edge(b,_)| (*a,*b))
                    })
                );
            }
        }
        deleted_edges.drain().for_each(|(a,b)| self.delete_edge(a,b));
    }

    pub fn delete_edge(&mut self, from: HaplotypeId, to:HaplotypeId) {
        self.succ.get_mut(&from).unwrap().retain(|Edge(succ,_)| succ != &to);
        self.pred.get_mut(&to).unwrap().retain(|Edge(pred,_)| pred != &from);
    }

    pub fn scaffold_haplotypes(&mut self, variant_positions:&HashSet<(usize,usize)>, min_reads:usize, min_frac:f64, min_snv:usize) -> HashMap<HaplotypeId,Haplotype> {
        assert!(min_frac > 0.5);
        self.remove_weak_edges(min_reads, min_frac);
        assert!(self.succ.keys().all(|hid| self.out_degree(hid) < 2 && self.in_degree(hid) < 2));
        let mut scaffolds: HashMap<HaplotypeId,Haplotype> = HashMap::new();
        for ht_from in self.succ.keys().filter(|&hid| self.in_degree(hid) == 0) {
            let mut haplotype_nodes = self.scaffold_from(*ht_from);
            assert!(!haplotype_nodes.is_empty());
            let ht_name = unsafe { haplotype_nodes.pop().unwrap_unchecked() };
            let mut scf = self.haplotypes[&ht_name].clone(); // TODO: do I really need to clone it?
            while let Some(ht_name) = haplotype_nodes.pop() {
                scf.extend(self.haplotypes[&ht_name].clone()); // TODO: same as above, can I safely remove it from the hashmap and not clone it?
            }
            scf.trim(variant_positions,0);
            assert!(scf.size() > 0);
            if scf.size() >= min_snv {
                scaffolds.insert(scf.uid(), scf);
            }
        }
        scaffolds
    }

    // Prerequisite: linear path from hid in the graph, no loops
    fn scaffold_from(&self, mut hid:HaplotypeId) -> Vec<HaplotypeId> {
        let mut haplotype_nodes = vec![];
        while self.out_degree(&hid) > 0 {
            haplotype_nodes.push(hid);
            hid = unsafe { self.succ[&hid].first().unwrap_unchecked().0 };
        }
        haplotype_nodes.push(hid);
        haplotype_nodes.reverse();
        haplotype_nodes
    }

    pub fn write_dot(&self, path:&Path) -> std::io::Result<()> {
        let file = File::create(path)?;
        
        let mut dot = BufWriter::new(file);
        dot.write_all(b"digraph \"\" {\n")?;
        dot.write_all(b"\tgraph [rankdir=LR, splines=true];\n")?;
        dot.write_all(b"\tnode [label=\"\\N\"];\n")?;

        for node in self.succ.keys() {
            let node_line = format!("\t\"{node}\" [style=filled, fillcolor=white];\n");
            dot.write_all(node_line.as_bytes())?;
        }

        for (from, adj) in self.succ.iter() {
            for Edge(to,count) in adj {
                let edge_line = format!("\t\"{from}\" -> \"{to}\" [label={count}];\n");
                dot.write_all(edge_line.as_bytes())?;
            }
        }

        dot.write_all(b"}\n")
    }
}
