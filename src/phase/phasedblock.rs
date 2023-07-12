use super::haplotree::HaploTree;

pub struct PhasedBlock {
    tid: usize,
    haplotree: HaploTree,
    haplotypes: Vec<usize>, // placeholder

}

impl PhasedBlock {
    pub fn new(tid:usize) -> PhasedBlock {
        PhasedBlock { 
            tid: tid,
            haplotree: HaploTree::new(),
            haplotypes: vec![]
        }
    }
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

//     def remove_haplotype(self, ht_id):
//         self.haplotypes.pop(ht_id,None)
//         self.haplotype_node.pop(ht_id,None)

//     def extend(self, snv_position, edges) -> bool:
//         is_ambiguous = False
//         self.end = snv_position
//         for ht_id in list(self.haplotypes.keys()): # TODO: use .items() ?
//             ht = self.haplotypes[ht_id]
//             ht_parent = self.haplotype_node[ht_id]
//             ht_nucleotide = ht.sequence[-1]
//             successors = (nuc_to for nuc_from,nuc_to in edges if nuc_from == ht_nucleotide)
//             snv = (snv_position,next(successors))
//             ht.append(snv)
//             self.haplotype_node[ht_id] = self.haplotree.extend(ht_parent,snv)
//             nuc_to = next(successors,None)
//             while nuc_to:
//                 is_ambiguous = True
//                 snv = (snv_position,nuc_to)
//                 new_node = self.haplotree.extend(ht_parent,snv)
//                 new_ht = copy.deepcopy(ht)
//                 new_ht.id = new_node
//                 new_ht.sequence = ht.sequence[:-1] + snv[1]
//                 self.haplotypes[new_node] = new_ht
//                 self.haplotype_node[new_node] = new_node
//                 nuc_to = next(successors,None)
//         return is_ambiguous

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