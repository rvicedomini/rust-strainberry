use itertools::Itertools;
use rayon::prelude::*;

use super::alignment::Alignment;
use super::window::Window;

const WINDOW_LEN: usize = 500;

pub struct Polisher<'a> {
    ref_sequences: &'a[Vec<u8>],
    read_sequences: &'a[Vec<u8>],
    alignments: Vec<Alignment>,
    windows: Vec<Window>,
    id_to_first_window_id: Vec<usize>,
    // targets_coverages: Vec<usize>
}

impl<'a> Polisher<'a> {

    pub fn new(ref_sequences: &'a [Vec<u8>], read_sequences: &'a [Vec<u8>], alignments: Vec<Alignment>) -> Self {
        Self {
            ref_sequences,
            read_sequences,
            alignments,
            windows: Vec::new(),
            id_to_first_window_id: Vec::new(),
            // targets_coverages: Vec::new(),
        }
    }

    pub fn initialize(&mut self) {
        
        self.find_alignment_breaking_points();
        self.initialize_windows();
        self.add_window_layers();
    }

    pub fn find_alignment_breaking_points(&mut self) {

        self.alignments.par_iter_mut().for_each(|a| {
            a.find_breaking_points(self.ref_sequences, self.read_sequences, WINDOW_LEN).unwrap();
        });
    }

    pub fn initialize_windows(&mut self) {

        self.windows.clear();
        // self.targets_coverages = vec![0; self.ref_sequences.len()];

        self.id_to_first_window_id = vec![0; self.ref_sequences.len()+1];
        for (id, seq) in self.ref_sequences.iter().enumerate() {
            let mut rank = 0;
            for j in (0..seq.len()).step_by(WINDOW_LEN) {
                let length = std::cmp::min(j + WINDOW_LEN, seq.len()) - j;
                self.windows.push(Window::build(id, rank, &seq[j..j+length]));
                rank += 1;
            }
            self.id_to_first_window_id[id+1] = self.id_to_first_window_id[id] + rank;
        }
    }

    pub fn add_window_layers(&mut self) {

        for a in std::mem::take(&mut self.alignments) {
            
            // self.targets_coverages[a.target_idx] += 1;

            let sequence = self.read_sequences[a.query_idx].as_slice();
            let breaking_points = &a.breaking_points;

            for (first, last) in breaking_points.iter().tuples() {
                if 100 * (last.1 - first.1) < 2 * WINDOW_LEN { // last.1 - first.1 < 0.02 * opts.window_len
                    continue
                }

                let window_id = self.id_to_first_window_id[a.target_idx] + first.0 / WINDOW_LEN;
                let window_start = (first.0 / WINDOW_LEN) * WINDOW_LEN;
                let window_seq = if a.strand == b'+' {
                    sequence[first.1..last.1].to_vec()
                } else {
                    let (qry_beg, qry_end) = (a.query_len-last.1, a.query_len-first.1);
                    crate::racon::sequence::revcomp(sequence[qry_beg..qry_end].to_vec())
                };

                self.windows[window_id].add_layer(window_seq, first.0 - window_start, last.0 - window_start - 1);
            }
        }
    }

    pub fn polish(&mut self) -> Vec<Vec<u8>> {

        self.windows.par_iter_mut().for_each(|win| {
            win.generate_consensus();
        });

        let mut result = Vec::new();
        let mut polished = Vec::new();
        for (i, win) in self.windows.iter().enumerate() {
            if let Some(consensus) = win.consensus.as_ref() {
                polished.extend(consensus);
            }
            if i == self.windows.len()-1 || self.windows[i+1].rank == 0 {
                result.push(std::mem::take(&mut polished));
            }
        }

        result
    }
    
}

