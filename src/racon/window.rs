#[derive(Debug, Default)]
pub struct Window {
    pub id: usize,
    pub rank: usize,
    pub consensus: Option<Vec<u8>>,
    pub sequences: Vec<Vec<u8>>,
    pub positions: Vec<(usize, usize)>,

}

impl Window {
    
    pub fn build(id:usize, rank:usize, backbone:&[u8]) -> Self {
        Self {
            id,
            rank,
            consensus: None,
            sequences: vec![backbone.to_vec()],
            positions: vec![(0,0)]
        }
    }

    pub fn add_layer(&mut self, sequence:Vec<u8>, begin:usize, end: usize) {
        
        if sequence.is_empty() || begin == end {
            return
        }

        if begin >= end || begin > self.sequences[0].len() || end > self.sequences[0].len() {
            panic!("layer begin and end positions are invalid!");
        }

        self.sequences.push(sequence);
        self.positions.push((begin,end));
    }

    pub fn generate_consensus(&mut self) {
        
        if self.sequences.len() < 3 {
            self.consensus = Some(self.sequences[0].to_vec());
            return
        }

        let mut spoa_engine = spoa_rs::AlignmentEngine::new_linear(spoa_rs::AlignmentType::kNW, 3, -5, -4);
        let mut graph = spoa_rs::Graph::new();

        let seq = std::str::from_utf8(&self.sequences[0]).unwrap();
        let aln = spoa_engine.align(seq, &graph).1;
        graph.add_alignment(aln, seq);

        
        // TODO: add subgraph and update_alignment bindings to my cutomized crate spoa_rs
        let window_depth = self.sequences.len()-1;
        let backbone_len = self.sequences[0].len();
        let offset = self.sequences[0].len() / 10;
        for (i, seq) in self.sequences.drain(..).enumerate().skip(1) {
            if self.positions[i].0 < offset && self.positions[i].1 > backbone_len - offset {
                let seq = unsafe { std::str::from_utf8_unchecked(&seq) };
                let aln = spoa_engine.align(seq, &graph).1;
                graph.add_alignment(aln, seq);
            }
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
                spdlog::warn!("contig {} might be chimeric in window {}!", self.id, self.rank);
            }
        }
        
        self.consensus = Some(consensus);
    }

}