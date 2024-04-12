/*use tinyvec::{tiny_vec,TinyVec};


fn flip_strand(strand:u8) -> u8 {
    assert!(strand == b'+' || strand == b'-');
    if strand == b'+' { b'-' } else { b'+' }
}


struct Node {
    id: usize,
    edges: TinyVec<[usize;10]>,
    successors: TinyVec<[usize;10]>,
    // leaf: bool,
    // complete: bool,
}


pub struct AwareGraph {
    nodes: Vec<Node>,
    edges: Vec<BiEdge>,
    transitives: Vec<BiEdge>,

}

    def __init__(self, seq_dict, read_dict):
        self.nodes : Dict[tuple,Node] = dict()
        self.next_node_id = 0
        self.edges : Dict[tuple,BiEdge] = dict()
        self.transitives : Dict[tuple,TransEdge] = dict()
        self.visited = set()
        self.seq_dict = seq_dict
        self.read_dict = read_dict




#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Snv {
    pub pos: usize,
    pub nuc: u8,
}



impl Node {
    fn new(snv:Option<Snv>, parent:usize) -> Self {
        Self { snv, parent, successors: tiny_vec!() }
    }
}

impl Default for Node {
    fn default() -> Self {
        Self::new(None, 0)
    }
}

pub struct HaploTree {
    root: usize,
    nodes: Vec<Node>,
}

impl HaploTree {

    pub fn new() -> HaploTree {
        HaploTree { root: 0, nodes: vec![Node::default()] }
    }

    pub fn get_root(&self) -> usize {
        self.root
    }

    pub fn clear(&mut self) {
        self.nodes.clear();
        self.nodes.push(Node::default());
    }

    pub fn get_parent(&self, node:usize) -> usize {
        self.nodes[node].parent 
    }

    // pub fn get_snv(&self, node:usize) -> Option<Snv> { self.nodes[node].snv }

    // pub fn is_leaf(&self, node:usize) -> bool { self.nodes[node].leaf }

    // pub fn is_complete(&self, node:usize) -> bool { self.nodes[node].complete }

    // pub fn set_leaf(&mut self, node:usize, value:bool) { self.nodes[node].leaf = value; }

    // pub fn set_complete(&mut self, node:usize, value:bool) { self.nodes[node].complete = value; }

    // pub fn peek(&self, node:Option<usize>, snv:Snv) -> Option<usize> {
    //     let succ = self.nodes[node?].successors.iter()
    //         .filter(|&&succ| self.nodes[succ].snv == Some(snv))
    //         .next()?;
    //     Some(*succ)
    // }

    pub fn extend(&mut self, node:usize, snv:Snv) -> usize {

        if let Some(&succ) = self.nodes[node].successors.iter().find(|&&succ| self.nodes[succ].snv == Some(snv)) {
            return succ;
        }

        let succ = self.nodes.len();
        self.nodes.push(Node::new(Some(snv), node));
        self.nodes[node].successors.push(succ);

        succ
    }

    pub fn path_extend(&mut self, node:usize, path:&[Snv]) -> usize {
        path.iter().fold(node, |node,&snv| self.extend(node,snv))
    }
}*/