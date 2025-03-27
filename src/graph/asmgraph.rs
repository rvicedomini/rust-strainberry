use ahash::AHashMap as HashMap;
use tinyvec::{tiny_vec, TinyVec};

use super::biedge::{self, EdgeKey};


#[derive(Debug)]
pub struct AsmEdge {
    pub key: EdgeKey, // TODO: check whether I really need to store the key in the struct
    // pub observations: usize,
    // pub gaps: Vec<i32>,
    // self.gapseq = {} # defaultdict(list)
}

impl AsmEdge {

    pub fn new(key:EdgeKey) -> Self {
        Self {
            key
        }
    }
}


#[derive(Debug)]
pub struct AsmNode {
    pub id: usize,
    pub sequence: String,
    pub edges: TinyVec<[EdgeKey;10]>,
}

impl AsmNode {

    pub fn new(id:usize, sequence:String) -> Self {
        Self {
            id,
            sequence,
            edges: tiny_vec![],
        }
    }
}


#[derive(Debug, Default)]
pub struct AsmGraph {
    nodes: HashMap<usize,AsmNode>,
    edges: HashMap<EdgeKey,usize>,
    edge_data: Vec<AsmEdge>,
    next_node_id: usize,
}

impl AsmGraph {
    
    pub fn new() -> Self {
        Self {
            nodes: HashMap::default(),
            edges: HashMap::default(),
            edge_data: vec![],
            next_node_id: 0
        }
    }

    pub fn nb_nodes(&self) -> usize { self.nodes.len() }
    pub fn nb_edges(&self) -> usize { self.edges.len() }

    pub fn len(&self) -> usize { self.nb_nodes() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    fn contains_edge(&self, edge_key: &EdgeKey) -> bool {
        let edge_key = biedge::canonical_edgekey(edge_key);
        self.edges.contains_key(&edge_key)
    }

    fn get_biedge_idx(&self, idx: usize) -> &AsmEdge { &self.edge_data[idx] }
    fn get_biedge_idx_mut(&mut self, idx: usize) -> &mut AsmEdge { &mut self.edge_data[idx] }

    fn get_biedge(&self, edge_key:&EdgeKey) -> Option<&AsmEdge> {
        let edge_key = biedge::canonical_edgekey(edge_key);
        let edge_idx = *self.edges.get(&edge_key)?;
        self.edge_data.get(edge_idx)
    }

    fn get_biedge_mut(&mut self, edge_key:&EdgeKey) -> Option<&mut AsmEdge> {
        let edge_key = biedge::canonical_edgekey(edge_key);
        let edge_idx = *self.edges.get(&edge_key)?;
        self.edge_data.get_mut(edge_idx)
    }

    fn get_biedge_or_create(&mut self, edge_key:&EdgeKey) -> &mut AsmEdge {
        let edge_key= biedge::canonical_edgekey(edge_key).into_owned();
        let edge_id = *self.edges.entry(edge_key).or_insert_with_key(|edge_key| {
            let node_from = self.nodes.get_mut(&edge_key.id_from).unwrap();
            node_from.edges.push(*edge_key);
            let node_to = self.nodes.get_mut(&edge_key.id_to).unwrap();
            node_to.edges.push(biedge::flip_edgekey(edge_key));
            let new_edge_id = self.edge_data.len();
            self.edge_data.push(AsmEdge::new(*edge_key));
            new_edge_id
        });
        unsafe { self.edge_data.get_mut(edge_id).unwrap_unchecked() }
    }

}


// class Node:
    
//     def __init__(self,seq_id, seq='', data=None) -> None:
//         self.seq_id = seq_id
//         self.sequence = seq
//         self.data = data
//         self.edges = dict()
//         self.transitives = dict()

//     def length(self) -> int:
//         return len(self.sequence)
    
//     def is_phased(self) -> bool:
//         return self.seq_id.endswith('_h')


// class BiEdge:
    
//     def __init__(self, biedge_key, data=None) -> None:
//         self.key = canon_key(biedge_key)
//         self.data = data

// class TransEdge:

//     def __init__(self, edge_key, data=None) -> None:
//         self.key = canon_key(edge_key)
//         self.data = data


// class AsmGraph(object):

//     def __init__(self):
//         self.nodes = dict()
//         self.edges = dict()
//         self.transitives = dict()
//         self.visited = set()

//     def add_node(self, seq_id, seq='', data=None):
//         self.nodes[seq_id] = Node(seq_id, seq, data)

//     def add_edge(self, key, data=None) -> BiEdge:
//         biedge = self._get_biedge_or_create(key)
//         biedge.data = data
//         return biedge
    
//     def get_edge(self, edge_key) -> BiEdge:
//         return self.edges[canon_key(edge_key)]

//     def _get_biedge_or_create(self, edge_key):
//         edge_key = canon_key(edge_key)
//         if edge_key in self.edges:
//             return self.edges[edge_key]
//         # create edge
//         edge = BiEdge(edge_key)
//         self.edges[edge_key] = edge
//         # update a and b nodes
//         a_id,a_dir,b_id,b_dir = edge_key
//         self.nodes[a_id].edges[edge_key] = edge
//         self.nodes[b_id].edges[(b_id,b_dir,a_id,a_dir)] = edge
//         return edge

//     def edge_exists(self, edge_key):
//         edge_key = canon_key(edge_key)
//         return (edge_key in self.edges)
    
//     def out_edges(self, node_id, node_dir) -> list:
//         return [key for key in self.nodes[node_id].edges.keys() if key[1]==node_dir]
    
//     def add_transitive(self, key, data=None) -> BiEdge:
//         edge = self._get_transitive_or_create(key)
//         edge.data = data
//         return edge
    
//     def get_transitive(self, edge_key) -> BiEdge:
//         return self.transitives[canon_key(edge_key)]

//     def _get_transitive_or_create(self, edge_key):
//         edge_key = canon_key(edge_key)
//         if edge_key in self.transitives:
//             return self.transitives[edge_key]
//         # create edge
//         edge = BiEdge(edge_key)
//         self.transitives[edge_key] = edge
//         # update a and b nodes
//         a_id,a_dir,b_id,b_dir = edge_key
//         self.nodes[a_id].transitives[edge_key] = edge
//         self.nodes[b_id].transitives[(b_id,b_dir,a_id,a_dir)] = edge
//         return edge

//     def transitive_exists(self, edge_key):
//         edge_key = canon_key(edge_key)
//         return (edge_key in self.transitives)
    
//     def out_transitives(self, node_id, node_dir) -> list:
//         return [key for key in self.nodes[node_id].transitives.keys() if key[1]==node_dir]

//     def remove_nodes_from(self, iterable):
//         for node_id in iterable:
//             self.nodes.pop(node_id,None)
    
//     def remove_edge(self, edge_key):
//         edge_key = canon_key(edge_key)
//         edge = self.edges[edge_key]
//         # remove edge pointers from nodes
//         self.nodes[edge_key[0]].edges.pop(edge_key,None)
//         flipped_key = flip_key(edge_key)
//         self.nodes[flipped_key[0]].edges.pop(flipped_key,None)
//         # remove edge from main dictionary
//         self.edges.pop(edge_key,None)

//     def remove_edges_from(self, iterable):
//         for edge_key in iterable:
//             self.remove_edge(edge_key)

//     def remove_transitive(self, edge_key):
//         edge_key = canon_key(edge_key)
//         edge = self.transitives[edge_key]
//         # remove edge pointers from nodes
//         self.nodes[edge_key[0]].transitives.pop(edge_key,None)
//         flipped_key = flip_key(edge_key)
//         self.nodes[flipped_key[0]].transitives.pop(flipped_key,None)
//         # remove edge from main dictionary
//         self.transitives.pop(edge_key,None)

//     def remove_transitives_from(self, iterable):
//         for edge_key in iterable:
//             self.remove_transitive(edge_key)

//     def generate_contigs(self):
//         # remove edges between unphased sequences (likely a SV in between)
//         # removed_edges = []
//         # for key,edge in self.edges.items():
//         #     from_deg = sum(1 for _,from_dir,_,_ in self.nodes[key[0]].edges.keys() if from_dir == key[1])
//         #     to_deg = sum(1 for _,to_dir,_,_ in self.nodes[key[2]].edges.keys() if to_dir == key[3])
//         #     if from_deg == 1 and to_deg == 1:
//         #         removed_edges.append(key)
//         # logger.debug(f'removing indel (?) edges: {removed_edges}')
//         # self.remove_edges_from(removed_edges)

//         contigs = []
//         visited = set()
//         # process sequences by decreasing length, starting from phased ones
//         for node in sorted(self.nodes.values(), key=lambda x:(x.is_phased(),x.length()), reverse=True):
//             if ((node.seq_id,'+') in visited) or ((node.seq_id,'-') in visited):
//                 continue
//             ctg_id = node.seq_id
//             ctg_sequence = ''
//             ctg_info = set()

//             in_edges = list(key for key in node.edges.keys() if key[1]=='+')
//             if len(in_edges) == 1 and (self.nodes[in_edges[0][0]].is_phased()) and (not self.nodes[in_edges[0][2]].is_phased()) and (in_edges[0][2],in_edges[0][3]) not in visited:
//                 visited.add((in_edges[0][2],flip_strand(in_edges[0][3])))
//                 ctg_sequence += (self.nodes[in_edges[0][2]].sequence if in_edges[0][3] == '-' else reverse_complement(self.nodes[in_edges[0][2]].sequence))
//                 ctg_info.update(self.nodes[in_edges[0][2]].data)
//                 edge = self.get_edge(in_edges[0])
//                 edge_gaplen = edge.data['gap']
//                 edge_gapseq = edge.data['seq']
//                 if edge_gaplen > 0 and (in_edges[0][2],in_edges[0][3]) in edge_gapseq:
//                     ctg_sequence += edge_gapseq[(in_edges[0][2],in_edges[0][3])]
//                 elif edge_gaplen > 0:
//                     ctg_sequence += ('N' * edge_gaplen)
//                 elif edge_gaplen < 0:
//                     ctg_sequence = ctg_sequence[-edge_gaplen:]
//             elif len(in_edges) == 1 and (self.nodes[in_edges[0][0]].is_phased()) and (not self.nodes[in_edges[0][2]].is_phased()):
//                 edge = self.get_edge(in_edges[0])
//                 edge_gaplen = edge.data['gap']
//                 edge_gapseq = edge.data['seq']
//                 if edge_gaplen > 0 and (in_edges[0][2],in_edges[0][3]) in edge_gapseq and (not edge_gapseq[(in_edges[0][2],in_edges[0][3])].startswith('N')):
//                     ctg_sequence += edge_gapseq[(in_edges[0][2],in_edges[0][3])]

//             ctg_sequence += node.sequence
//             ctg_info.update(node.data)

//             out_edges = list(key for key in node.edges.keys() if key[1]=='-')
//             if len(out_edges) == 1 and (self.nodes[out_edges[0][0]].is_phased()) and (not self.nodes[out_edges[0][2]].is_phased()) and (out_edges[0][2],out_edges[0][3]) not in visited:
//                 visited.add((out_edges[0][2],flip_strand(out_edges[0][3])))
//                 edge = self.get_edge(out_edges[0])
//                 edge_gaplen = edge.data['gap']
//                 edge_gapseq = edge.data['seq']
//                 if edge_gaplen > 0 and (out_edges[0][0],out_edges[0][1]) in edge_gapseq:
//                     ctg_sequence += edge_gapseq[(out_edges[0][0],out_edges[0][1])]
//                 elif edge_gaplen > 0:
//                     ctg_sequence += ('N' * edge_gaplen)
//                 elif edge_gaplen < 0:
//                     ctg_sequence = ctg_sequence[:edge_gaplen]
//                 ctg_sequence += self.nodes[out_edges[0][2]].sequence if out_edges[0][3] == '+' else reverse_complement(self.nodes[out_edges[0][2]].sequence)
//                 ctg_info.update(self.nodes[out_edges[0][2]].data)
//             elif len(out_edges) == 1 and (self.nodes[out_edges[0][0]].is_phased()) and (not self.nodes[out_edges[0][2]].is_phased()):
//                 edge = self.get_edge(out_edges[0])
//                 edge_gaplen = edge.data['gap']
//                 edge_gapseq = edge.data['seq']
//                 if edge_gaplen > 0 and (out_edges[0][0],out_edges[0][1]) in edge_gapseq and (not edge_gapseq[(out_edges[0][0],out_edges[0][1])].startswith('N')):
//                     ctg_sequence += edge_gapseq[(out_edges[0][0],out_edges[0][1])]
            
//             contigs.append((ctg_id,ctg_sequence,ctg_info))
//         return contigs


//     def load_from_gfa(self, filename, load_seq=True):
//         segments = []
//         links = []
//         with open(filename,'r') as gfa:
//             for line in gfa:
//                 if line.startswith('S'):
//                     fields = line.strip().split('\t')
//                     name = fields[1]
//                     sequence = fields[2] if (load_seq and fields[2] != '*') else ''
//                     segments.append((name,sequence))
//                 elif line.startswith('L'):
//                     fields = line.strip().split('\t')
//                     from_name,from_dir,to_name,to_dir,overlap = fields[1:6]
//                     links.append((from_name,flip_strand(from_dir),to_name,to_dir))
//         for name,sequence in segments:
//             self.add_node(name,sequence)
//         for key in links:
//             if (key[0] in self.nodes) and (key[2] in self.nodes):
//                 self.add_edge(key)

//     def write_gfa(self, gfa_path, write_seq=False):
//         with open(f'{gfa_path}','w') as gfa:
//             gfa.write(f'H\tVN:Z:1.0\n')
//             for node_id in self.nodes.keys():
//                 node_seq = self.nodes[node_id].sequence if (write_seq and len(self.nodes[node_id].sequence) > 0) else '*'
//                 node_length = len(self.nodes[node_id].sequence)
//                 gfa.write(f'S\t{node_id}\t{node_seq}\tLN:i:{node_length}\n')
//             for edge in self.edges.values():
//                 first_id,first_dir,second_id,second_dir = edge.key
//                 gfa.write(f'L\t{first_id}\t{flip_strand(first_dir)}\t{second_id}\t{second_dir}\t0M\n')

//     def write_dot(self, path):
//         with open(f'{path}','w') as dot:
//             dot.write('digraph "" {\n')
//             dot.write('\tgraph [rankdir=LR, splines=true];\n')
//             # dot.write('\tnode [label="\\N"];\n')
//             for node_id in self.nodes.keys():
//                 fillcolor = 'orange' if node_id.endswith('_h') else 'white'
//                 dot.write(f'\t{node_id}\t[label="{node_id}", style=filled, fillcolor={fillcolor}];\n')
//             for edge in self.edges.values():
//                 first_id,first_dir,second_id,second_dir = edge.key
//                 arrowtail = 'normal' if first_dir=='+' else 'inv'
//                 arrowhead = 'normal' if second_dir=='+' else 'inv'
//                 gap_len = edge.data['gap'] if (edge.data and 'gap' in edge.data) else ''
//                 dot.write(f'\t{first_id} -> {second_id}\t[arrowtail={arrowtail}, arrowhead={arrowhead}, dir=both, label="{gap_len}"];\n')
//             for edge in self.transitives.values():
//                 first_id,first_dir,second_id,second_dir = edge.key
//                 arrowtail = 'normal' if first_dir=='+' else 'inv'
//                 arrowhead = 'normal' if second_dir=='+' else 'inv'
//                 data = edge.data
//                 dot.write(f'\t{first_id} -> {second_id}\t[arrowtail={arrowtail}, arrowhead={arrowhead}, dir=both, color="red", label="{data:.4f}"];\n')
//             dot.write('}')

//     def write_fasta(self, fasta_path):
//         with open(fasta_path,'w') as fasta:
//             for node in self.nodes.values():
//                 fasta.write(f'>{node.seq_id}\n{insert_newlines(node.sequence)}\n')
