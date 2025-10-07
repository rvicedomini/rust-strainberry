use itertools::Itertools;

use super::biedge::EdgeKey;

pub struct Junction {
    pub in_edges: Vec<EdgeKey>,
    pub mid_edges: Vec<EdgeKey>,
    pub out_edges: Vec<EdgeKey>,
}

impl Junction {

    pub fn new() -> Self {
        Self {
            in_edges: vec![],
            mid_edges: vec![],
            out_edges: vec![],
        }
    }

    pub fn size(&self) -> usize {
        self.mid_edges.len() + 1
    }

    pub fn input_size(&self) -> usize {
        self.in_edges.len()
    }

    pub fn output_size(&self) -> usize {
        self.out_edges.len()
    }

    pub fn input_nodes(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.in_edges.iter().map(|key| key.dst_id())
    }

    pub fn inner_nodes(&self) -> impl Iterator<Item = usize> + use<'_> {
        std::iter::once(self.in_edges[0].src_id())
            .chain(self.mid_edges.iter().map(|key| key.dst_id()))
    }

    pub fn output_nodes(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.out_edges.iter().map(|key| key.dst_id())
    }

    pub fn inout_nodes(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.input_nodes().chain(self.output_nodes())
    }

    pub fn input_dir(&self) -> u8 {
        self.in_edges[0].src_dir()
    }

    pub fn output_dir(&self) -> u8 {
        self.out_edges[0].src_dir()
    }

    pub fn inputs(&self) -> impl Iterator<Item = (usize,u8)> + use<'_> {
        self.in_edges.iter().map(|key| key.dst())
    }

    pub fn outputs(&self) -> impl Iterator<Item = (usize,u8)> + use<'_> {
        self.out_edges.iter().map(|key| key.dst())
    }

}


impl std::fmt::Display for Junction {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let in_nodes = self.input_nodes().join(",");
        let mid_nodes = self.inner_nodes().join(",");
        let out_nodes = self.output_nodes().join(",");
        write!(f, "Junction({in_nodes}|{mid_nodes}|{out_nodes})")
    }
}
