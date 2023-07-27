use tinyvec::{tiny_vec,TinyVec};

#[derive(Clone,Copy)]
pub struct SNV {
    pub pos: usize,
    pub nuc: u8,
}

struct Node {
    snv: Option<SNV>,
    parent: usize,
    successors: TinyVec<[usize;10]>,
    leaf: bool,
    complete: bool,
}

impl Node {
    fn default() -> Node {
        Node {
            snv: None,
            parent: 0,
            successors: tiny_vec!(),
            leaf: false,
            complete: false
        }
    }
}

pub struct HaploTree {
    root: usize,
    nodes: Vec<Node>,
}

impl HaploTree {

    pub fn new() -> HaploTree {
        let root_node = Node::default();
        HaploTree { root: 0, nodes: vec![root_node] }
    }

    pub fn get_root(&self) -> usize { self.root }

    pub fn clear(&mut self) {
        self.nodes.clear();
        self.nodes.push(Node::default());
    }

    pub fn get_parent(&self, node:usize) -> usize { self.nodes[node].parent }

    pub fn get_snv(&self, node:usize) -> Option<SNV> { self.nodes[node].snv }

    pub fn is_leaf(&self, node:usize) -> bool { self.nodes[node].leaf }

    pub fn is_complete(&self, node:usize) -> bool { self.nodes[node].complete }

    pub fn set_leaf(&mut self, node:usize, value:bool) { self.nodes[node].leaf = value; }

    pub fn set_complete(&mut self, node:usize, value:bool) { self.nodes[node].complete = value; }

    pub fn peek(&self, node:Option<usize>, snv:SNV) -> Option<usize> {

        if let Some(node) = node {
            if let Some(succ) = self.nodes[node].successors.iter().filter(|&succ| matches!(self.nodes[*succ].snv,snv)).next() {
                return Some(*succ)
            }
        }
        None
    }

    pub fn extend(&mut self, node:usize, snv:SNV) -> usize {

        if let Some(succ) = self.nodes[node].successors.iter().filter(|&succ| matches!(self.nodes[*succ].snv,snv)).next() {
            return *succ
        }

        let succ = self.nodes.len();
        self.nodes.push(Node {
            snv: Some(snv),
            parent: node,
            successors: tiny_vec!(),
            leaf: false,
            complete: false 
        });
        
        let node_successors = &mut self.nodes[node].successors;
        if !node_successors.contains(&succ) { // I do not expect node_successor to be big
            node_successors.push(succ);
        }

        succ
    }

    pub fn path_extend(&mut self, node:usize, path:Vec<SNV>) -> usize {
        path.iter().fold(node, |node,&snv| self.extend(node,snv))
    }

    // pub fn write_dot(&self, path:Path) { todo!() }
}