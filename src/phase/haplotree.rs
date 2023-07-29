use tinyvec::{tiny_vec,TinyVec};

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct SNV {
    pub pos: usize,
    pub nuc: u8,
}

struct Node {
    snv: Option<SNV>,
    parent: usize,
    successors: TinyVec<[usize;10]>,
    // leaf: bool,
    // complete: bool,
}

impl Node {
    fn new(snv:Option<SNV>, parent:usize) -> Self {
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

    // pub fn get_snv(&self, node:usize) -> Option<SNV> { self.nodes[node].snv }

    // pub fn is_leaf(&self, node:usize) -> bool { self.nodes[node].leaf }

    // pub fn is_complete(&self, node:usize) -> bool { self.nodes[node].complete }

    // pub fn set_leaf(&mut self, node:usize, value:bool) { self.nodes[node].leaf = value; }

    // pub fn set_complete(&mut self, node:usize, value:bool) { self.nodes[node].complete = value; }

    // pub fn peek(&self, node:Option<usize>, snv:SNV) -> Option<usize> {
    //     let succ = self.nodes[node?].successors.iter()
    //         .filter(|&&succ| self.nodes[succ].snv == Some(snv))
    //         .next()?;
    //     Some(*succ)
    // }

    pub fn extend(&mut self, node:usize, snv:SNV) -> usize {

        if let Some(&succ) = self.nodes[node].successors.iter().filter(|&&succ| self.nodes[succ].snv == Some(snv)).next() {
            return succ;
        }

        let succ = self.nodes.len();
        self.nodes.push(Node::new(Some(snv), node));
        assert!(!self.nodes[node].successors.contains(&succ));
        self.nodes[node].successors.push(succ);

        succ
    }

    pub fn path_extend(&mut self, node:usize, path:&Vec<SNV>) -> usize {
        path.iter().fold(node, |node,&snv| self.extend(node,snv))
    }
}