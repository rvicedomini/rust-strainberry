use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use itertools::Itertools;
use rustc_hash::FxHashMap;

use super::haplotype::Haplotype;

type HaplotypeId = (usize,usize,usize,usize);

pub struct HaploGraph {
    graph: FxHashMap<HaplotypeId,Vec<(HaplotypeId,usize)>>
}

impl HaploGraph {

    pub fn new(haplotypes: &FxHashMap<HaplotypeId,Haplotype>, sread_haplotypes: &FxHashMap<String, Vec<HaplotypeId>>) -> HaploGraph {
        let mut graph = FxHashMap::default();
        
        for ht_uid in haplotypes.keys() {
            graph.insert(*ht_uid, vec![]);
        }

        for ht_list in sread_haplotypes.values() {
            for (a,b) in ht_list.iter().sorted_unstable_by_key(|&hid| haplotypes[hid].beg()).tuple_windows::<(_,_)>() {
                let a_adj = graph.get_mut(a).unwrap();
                if let Some(idx) = a_adj.iter().position(|(uid,_)| uid == b) {
                    let p = &mut a_adj[idx];
                    p.1 += 1;
                } else {
                    a_adj.push((*b,1));
                }
            }
        }

        HaploGraph{graph}
    }

    pub fn write_dot(&self, path:&Path) -> std::io::Result<()> {

        let file = File::create(path)?;
        let mut dot = BufWriter::new(file);

        dot.write_all(b"digraph \"\" {\n")?;
        dot.write_all(b"\tgraph [rankdir=LR, splines=true];\n")?;

        dot.write_all(b"\tnode [label=\"\\N\"];\n")?;

        for node in self.graph.keys() {
            let node_line = format!("\t\"{node:?}\" [style=filled, fillcolor=white];\n");
            dot.write_all(node_line.as_bytes())?;
        }

        for (a, succ) in self.graph.iter() {
            let a_name = format!("{a:?}");
            for (b,cnt) in succ.iter() {
                let b_name = format!("{b:?}");
                let edge_line = format!("\t\"{a_name}\" -> \"{b_name}\" [label={cnt}];\n");
                dot.write_all(edge_line.as_bytes())?;
            }
        }

        dot.write_all(b"}\n")?;

        Ok(())
    }
}


// def write_dot(self, path, id_as_name=True, output_transitives=True):
//     with open(f'{path}','w') as dot:
//         dot.write('digraph "" {\n')
//         dot.write('\tgraph [rankdir=LR, splines=true];\n')
//         # dot.write('\tnode [label="\\N"];\n')
//         for node_id, node in self.nodes.items():
//             node_name = f'{node_id}' if id_as_name else f'{node.ctg.name()}'
//             node_length = node.ctg.length()
//             node_depth = node.ctg.alignedbases/node_length
//             fillcolor = 'orange' if node.ctg.phased else 'white'
//             dot.write(f'\t{node_id}\t[label="{node_name} ({int(node_depth)}X)", style=filled, fillcolor={fillcolor}];\n')
//         for edge in self.edges.values():
//             first_id,first_dir,second_id,second_dir = edge.key
//             arrowtail = 'normal' if first_dir=='+' else 'inv'
//             arrowhead = 'normal' if second_dir=='+' else 'inv'
//             edge_gap = 0 if len(edge.gaps) == 0 else int(statistics.median(edge.gaps))
//             dot.write(f'\t{first_id} -> {second_id}\t[arrowtail={arrowtail}, arrowhead={arrowhead}, dir=both, label="{edge.observations}/{edge_gap}bp"];\n')
//         if output_transitives:
//             for tedge in self.transitives.values():
//                 first_id,first_dir,second_id,second_dir = tedge.key
//                 arrowtail = 'normal' if first_dir=='+' else 'inv'
//                 arrowhead = 'normal' if second_dir=='+' else 'inv'
//                 tedge_gap = 0 if len(tedge.gaps) == 0 else int(statistics.median(tedge.gaps))
//                 dot.write(f'\t{first_id} -> {second_id}\t[arrowtail={arrowtail}, arrowhead={arrowhead}, dir=both, color="red", label="{tedge.observations}/{tedge_gap}bp"];\n')
//         dot.write('}')