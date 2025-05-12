// use std::borrow::Cow;

use crate::bitseq::BitSeq;


#[derive(Debug)]
pub struct Layer<'a> {
    pub sequence: (&'a BitSeq, usize, usize),
    pub strand: u8,
    pub begin: usize,
    pub end: usize,
}

impl Layer<'_> {
    pub fn len(&self) -> usize { self.sequence.2 - self.sequence.1 }
    pub fn is_empty(&self) -> bool { self.len() == 0 }
}


#[derive(Debug, Default)]
pub struct Window <'a> {
    pub id: usize,
    pub rank: usize,
    pub consensus: Option<Vec<u8>>,
    pub layers: Vec<Layer<'a>>
}

impl<'a> Window<'a> {
    
    pub fn build(id:usize, rank:usize, backbone: (&'a BitSeq, usize, usize)) -> Self {
        Self {
            id,
            rank,
            consensus: None,
            layers: vec![ Layer { sequence:backbone, strand:b'+', begin:0, end:0 } ]
        }
    }

    pub fn add_layer(&mut self, sequence: (&'a BitSeq, usize, usize), strand: u8, begin: usize, end: usize) {
        
        if sequence.1 == sequence.2 || begin == end {
            return
        }

        if begin >= end || begin > self.layers[0].len() || end > self.layers[0].len() {
            panic!("layer begin and end positions are invalid!");
        }

        self.layers.push(Layer { sequence, strand, begin, end });
    }

    pub fn generate_consensus(&mut self) {

        let (sequence, beg, end) = self.layers[0].sequence;
        
        if self.layers.len() < 3 {
            self.consensus = Some(sequence.subseq(beg, end));
            return
        }

        let mut spoa_engine = spoa_rs::AlignmentEngine::new_linear(
            spoa_rs::AlignmentType::kNW, 3, -5, -4
        );
        
        let mut graph = spoa_rs::Graph::new();
        let seq = sequence.subseq(beg, end);
        let seq = std::str::from_utf8(&seq).unwrap();
        let alignment = spoa_engine.align(&graph, seq);
        graph.add_alignment(alignment, seq);

        self.layers[1..].sort_unstable_by_key(|l| l.begin);

        let window_depth = self.layers.len()-1;
        let backbone_len = self.layers[0].len();
        let offset = self.layers[0].len() / 100;
        for layer in self.layers.drain(..).skip(1) {

            let (sequence, beg, end) = layer.sequence;
            let mut seq = sequence.subseq(beg, end);
            if layer.strand == b'-' { crate::seq::revcomp_inplace(&mut seq); }

            // let seq = if layer.strand == b'+' {
            //     Cow::Borrowed(layer.sequence)
            // } else {
            //     Cow::Owned(crate::seq::revcomp(layer.sequence))
            // };
            
            let seq = unsafe { std::str::from_utf8_unchecked(&seq) };
            let alignment = if layer.begin < offset && layer.end > backbone_len - offset {
                spoa_engine.align(&graph, seq)
            } else {
                spoa_engine.align_subgraph(&graph, seq, layer.begin as u32, layer.end as u32)
            };
            graph.add_alignment(alignment, seq);
        }

        let (mut consensus, coverages) = graph.generate_consensus_with_coverage();
        
        let average_coverage = window_depth as u32 / 2;
        let begin = (0..coverages.len()).find(|i| coverages[*i] >= average_coverage);
        let end = (0..coverages.len()).rev().find(|i| coverages[*i] >= average_coverage);
        
        match (begin,end) {
            (Some(begin),Some(end)) if begin < end => {
                consensus = std::mem::take(&mut consensus[begin..end+1].to_vec());
            },
            _ => {
                spdlog::debug!("contig {} might be chimeric in window {}!", self.id, self.rank);
            }
        }
        
        self.consensus = Some(consensus);
    }

}