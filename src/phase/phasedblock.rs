use itertools::Itertools;
use rustc_hash::FxHashMap;

use crate::phase;

use super::haplotree::{HaploTree,SNV};
use super::haplotype::Haplotype;

pub struct PhasedBlock {
    tid: usize,
    haplotree: HaploTree,
    haplotypes: FxHashMap<usize,Haplotype>,
    haplotype_node: FxHashMap<usize,usize>,
    begin: usize,
}

impl PhasedBlock {

    pub fn new(tid:usize) -> PhasedBlock {
        PhasedBlock { 
            tid: tid,
            haplotree: HaploTree::new(),
            haplotypes: FxHashMap::default(),
            haplotype_node: FxHashMap::default(),
            begin: 0,
        }
    }

    pub fn init(&mut self, pos:usize, nucleotides:Vec<u8>) {
        
        self.haplotree.clear();
        self.haplotypes.clear();
        self.haplotype_node.clear();
        self.begin = pos;
        
        let htree_root = self.haplotree.get_root();
        for nuc in nucleotides {
            let snv = SNV{pos,nuc};
            let htree_node = self.haplotree.extend(htree_root, snv);
            self.haplotypes.insert(htree_node, Haplotype::new(htree_node, self.tid, vec![snv]));
            self.haplotype_node.insert(htree_node, htree_node);
        }
    }

    pub fn haplotypes(&self) -> &FxHashMap<usize,Haplotype> { &self.haplotypes }
    pub fn haplotypes_mut(&mut self) -> &mut FxHashMap<usize,Haplotype> { &mut self.haplotypes }

    pub fn get(&self, hid:usize) -> &Haplotype { &self.haplotypes[&hid] }

    pub fn begin(&self) -> usize { self.begin }

    pub fn remove_haplotype(&mut self, hid:usize) {
        self.haplotypes.remove(&hid);
        self.haplotype_node.remove(&hid);
    }

    pub fn drain(&mut self) -> Vec<Haplotype> {
        let mut haplotypes = Vec::with_capacity(self.haplotypes.len());
        for (i,mut ht) in self.haplotypes_mut().drain().map(|(_,ht)| ht).enumerate() {
            ht.set_hid(i);
            haplotypes.push(ht);
        }
        haplotypes
    }

    pub fn extend(&mut self, pos:usize, edges:Vec<(u8,u8)>) -> bool {
        let mut is_ambiguous = false;
        // self.end = pos;
        for hid in self.haplotypes.keys().cloned().collect_vec() {
            let ht = self.haplotypes.get_mut(&hid).unwrap();
            let ht_parent = self.haplotype_node[&hid];
            let ht_nuc = ht.last_nuc();
            let mut successors = edges.iter().filter(|(s,_)| *s == ht_nuc).map(|(_,t)| t);
            let snv = SNV{ pos, nuc: *successors.next().unwrap() };
            ht.push(snv);
            self.haplotype_node.insert(hid, self.haplotree.extend(ht_parent,snv));
            let mut new_haplotypes = vec![];
            while let Some(&nuc) = successors.next() {
                is_ambiguous = true;
                let snv = SNV{pos,nuc};
                let new_node = self.haplotree.extend(ht_parent,snv);
                let mut new_ht = ht.clone();
                new_ht.set_hid(new_node);
                new_ht.last_mut().nuc = nuc;
                new_haplotypes.push((new_node,new_ht));
            }
            for (new_node, new_ht) in new_haplotypes {
                self.haplotypes.insert(new_node,new_ht);
                self.haplotype_node.insert(new_node, new_node);
            }
        }
        is_ambiguous
    }

    pub fn split_and_init(&mut self, pos:usize, lookback:Option<usize>) -> PhasedBlock {

        let mut out_edges: FxHashMap<usize,Vec<usize>> = FxHashMap::default();
        for hid in self.haplotypes.keys() {
            let key = self.haplotree.get_parent(self.haplotype_node[hid]);
            out_edges.entry(key).or_default().push(*hid);
        }

        let mut phasedblock = PhasedBlock::new(self.tid);
        let new_root = phasedblock.haplotree.get_root();

        for mut ht_ids in out_edges.into_values() {
            while let Some(hid) = ht_ids.pop() {
                let ht = self.haplotypes.get_mut(&hid).unwrap();
                let mut new_ht = ht.split_off(ht.raw_size()-1, 0, lookback.unwrap_or(0));
                let new_hid = phasedblock.haplotree.path_extend(new_root, new_ht.raw_variants());
                new_ht.set_hid(new_hid);
                phasedblock.begin = new_ht.first_pos();
                phasedblock.haplotype_node.insert(new_hid, new_hid);
                phasedblock.haplotypes.insert(new_hid, new_ht);
                if ht_ids.len() > 0 {
                    self.haplotypes.remove(&hid);
                }
            }
        }

        phasedblock
    }

//     def split_and_init(self, snv_position, lookback=0) -> PhasedBlock:
//         out_edges = defaultdict(list)
//         for ht_id in self.haplotypes:
//             key = self.haplotree.get_parent(self.haplotype_node[ht_id])
//             out_edges[key].append(ht_id)
//         phased_block = PhasedBlock(self.reference_id, snv_position)
//         new_root = phased_block.haplotree.get_root()
        
//         for ht_list in out_edges.values():
//             while len(ht_list) > 0:
//                 ht_id = ht_list.pop()
//                 ht = self.haplotypes[ht_id]
//                 new_ht = ht.split_at(len(ht.positions)-1, 0, lookback)
//                 new_ht.id = phased_block.haplotree.path_extend(new_root,zip(new_ht.positions,new_ht.sequence))
//                 phased_block.haplotypes[new_ht.id] = new_ht
//                 phased_block.haplotype_node[new_ht.id] = new_ht.id
//                 phased_block.hap_begin = new_ht.positions[0]
//                 if len(ht_list) != 0:
//                     self.haplotypes.pop(ht_id,None)

}






// class PhasedBlock:
    
//     def __init__(self, reference_id, begin=None) -> None:
//         self.reference_id = reference_id
//         self.haplotree = HaploTree()
//         self.haplotypes = dict()
//         self.haplotype_node = dict()
//         self.hap_begin = begin
//         self.begin = begin
//         self.end = begin

//     def init_haplotypes(self, snv_position, nucleotides) -> None:
//         self.haplotree.clear()
//         self.haplotypes.clear()
//         self.haplotype_node.clear()
//         self.hap_begin = snv_position
//         self.begin = snv_position
//         self.end = snv_position
//         for nuc in nucleotides:
//             snv = (snv_position,nuc)
//             htree_root = self.haplotree.get_root()
//             htree_node = self.haplotree.extend(htree_root,snv)
//             self.haplotypes[htree_node] = NextHaplotype(htree_node,self.reference_id,[snv])
//             self.haplotype_node[htree_node] = htree_node

//     def split_and_init(self, snv_position, lookback=0) -> PhasedBlock:
//         out_edges = defaultdict(list)
//         for ht_id in self.haplotypes:
//             key = self.haplotree.get_parent(self.haplotype_node[ht_id])
//             out_edges[key].append(ht_id)
//         phased_block = PhasedBlock(self.reference_id, snv_position)
//         new_root = phased_block.haplotree.get_root()
        
//         for ht_list in out_edges.values():
//             while len(ht_list) > 0:
//                 ht_id = ht_list.pop()
//                 ht = self.haplotypes[ht_id]
//                 new_ht = ht.split_at(len(ht.positions)-1, 0, lookback)
//                 new_ht.id = phased_block.haplotree.path_extend(new_root,zip(new_ht.positions,new_ht.sequence))
//                 phased_block.haplotypes[new_ht.id] = new_ht
//                 phased_block.haplotype_node[new_ht.id] = new_ht.id
//                 phased_block.hap_begin = new_ht.positions[0]
//                 if len(ht_list) != 0:
//                     self.haplotypes.pop(ht_id,None)
                    
//         return phased_block