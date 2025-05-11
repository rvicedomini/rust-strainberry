use std::borrow::Cow;
use std::cmp::{PartialOrd, Eq};
use std::fmt::Display;
use std::hash::Hash;
use std::io::BufRead;
use std::path::Path;

use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use anyhow::{bail, Result};
use itertools::Itertools;
use tinyvec::{tiny_vec, TinyVec};

use crate::bitseq::BitSeq;
use crate::seq::flip_strand;


#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Link<NodeId> 
where
    NodeId: Hash + Eq + PartialOrd + Display + Clone
{
    pub id_from: NodeId,
    pub strand_from: u8,
    pub id_to: NodeId,
    pub strand_to: u8,
}

impl<NodeId> Link<NodeId>
where
    NodeId: Hash + Eq + PartialOrd + Display + Clone
{

    pub fn new(id_from: NodeId, strand_from: u8, id_to: NodeId, strand_to: u8) -> Self {
        Self { id_from, strand_from, id_to, strand_to }
    }

    pub fn flip(&mut self) {
        std::mem::swap(&mut self.id_from, &mut self.id_to);
        std::mem::swap(&mut self.strand_from, &mut self.strand_to);
    }

    pub fn is_canonical(&self) -> bool {
        self.id_from < self.id_to || (self.id_from == self.id_to && self.strand_from <= self.strand_to)
    }

    // returns whether the key was already canonical
    pub fn canonicalize(&mut self) -> bool {
        if !self.is_canonical() {
            self.flip();
            return false
        }
        true
    }

    #[inline(always)]
    pub fn canonical_link(link: &Link<NodeId>) -> Cow<'_, Link<NodeId>> {
        if link.is_canonical() {
            return Cow::Borrowed(link)
        }
        let mut link = link.clone();
        link.flip();
        Cow::Owned(link)
    }

}


#[derive(Debug)]
pub struct Segment<NodeId>
where
    NodeId: Hash + Eq + PartialOrd + Default + Display + Clone
{
    pub sequence: BitSeq,
    pub links: TinyVec<[Link<NodeId>;10]>
}

impl<NodeId> Segment<NodeId>
where
    NodeId: Hash + Eq + PartialOrd + Default + Display + Clone
{

    pub fn new(sequence:BitSeq) -> Self {
        Self {
            sequence,
            links: tiny_vec![],
        }
    }
}


#[derive(Debug, Default)]
pub struct AsmGraph<NodeId>
where
    NodeId: Hash + Eq + PartialOrd + Default + Display + Clone
{
    segments: HashMap<NodeId,Segment<NodeId>>,
    links: HashSet<Link<NodeId>>,
}

impl<NodeId> AsmGraph<NodeId>
where
    NodeId: Hash + Eq + PartialOrd + Default + Display + Clone
{
    
    pub fn new() -> Self {
        Self {
            segments: HashMap::new(),
            links: HashSet::new(),
        }
    }

    pub fn from_gfa(gfa_path: &Path) -> Result<AsmGraph<String>> {
        
        let mut graph: AsmGraph<String> = AsmGraph::new();

        let mut segments = Vec::new();
        let mut links = Vec::new();

        let gfa_reader = crate::utils::get_file_reader(gfa_path);
        for (line_num, line) in gfa_reader.lines().map_while(Result::ok).enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue
            }

            let fields = line.split('\t').collect_vec();
            match fields[0] {
                "S" if fields.len() >= 3 => {
                    let name = fields[1].to_string();
                    let sequence = fields[2].as_bytes().to_vec();
                    segments.push((name,sequence));
                },
                "L" if fields.len() >= 5 => {
                    let id_from = fields[1].to_string();
                    let strand_from = flip_strand(fields[2].bytes().next().unwrap());
                    let id_to = fields[3].to_string();
                    let strand_to = fields[4].bytes().next().unwrap();
                    links.push(Link::new(id_from, strand_from, id_to, strand_to));
                }
                "S"|"L" => { bail!("Cannot parse GFA at line {line_num}") }
                _ => {}
            }
        }

        segments.into_iter()
            .for_each(|(id,seg)| { graph.add_node(id.to_string(), BitSeq::from_utf8(&seg)); });

        links.into_iter()
            .for_each(|link| { graph.add_link(link); });

        Ok(graph)
    }

    pub fn nb_segments(&self) -> usize { self.segments.len() }
    pub fn nb_links(&self) -> usize { self.links.len() }

    pub fn len(&self) -> usize { self.nb_segments() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn add_node(&mut self, id:NodeId, sequence:BitSeq) {
        self.segments.insert(id, Segment::new(sequence));
    }

    pub fn contains_segment(&self, id: &NodeId) -> bool {
        self.segments.contains_key(id)
    }

    pub fn get_segment(&self, id: &NodeId) -> Option<&Segment<NodeId>> {
        self.segments.get(id)
    }

    pub fn segments(&self) -> std::collections::hash_map::Iter<'_, NodeId, Segment<NodeId>> {
        self.segments.iter()
    }

    pub fn contains_link(&self, link: &Link<NodeId>) -> bool {
        let canon_link = Link::canonical_link(link);
        self.links.contains(&canon_link)
    }

    pub fn add_link(&mut self, link: Link<NodeId>) -> bool {
        assert!(self.contains_segment(&link.id_from) && self.contains_segment(&link.id_to));
        let canon_link = Link::canonical_link(&link);
        self.links.insert(canon_link.into_owned())
    }

    /* OUTPUT METHODS */

    pub fn write_gfa(&self, gfa_path: &Path) -> std::io::Result<()> {

        let mut gfa = crate::utils::get_file_writer(gfa_path);
        gfa.write_all(b"H\tVN:Z:1.0\n")?;
        for (id, segment) in self.segments() {
            let seg_len = segment.sequence.len();
            let seg_line = if seg_len > 0 {
                format!("S\t{}\t{}\tLN:i:{}\n", id, segment.sequence, seg_len)
            } else {
                format!("S\t{}\t*\tLN:i:{}\n", id, seg_len)
            };
            gfa.write_all(seg_line.as_bytes())?;
        }

        for link in self.links.iter() {
            let Link { id_from, strand_from, id_to, strand_to } = link;
            let link_line = format!("L\t{}\t{}\t{}\t{}\t0M\n", id_from, flip_strand(*strand_from) as char, id_to, *strand_to as char);
            gfa.write_all(link_line.as_bytes())?;
        }

        Ok(())
    }

    pub fn write_fasta(&self, fasta_path: &Path, min_length:usize) -> std::io::Result<()> {

        let mut fasta_writer = crate::utils::get_file_writer(fasta_path);
        for (node_id, node) in self.segments() {
            if node.sequence.len() < min_length {
                continue;
            }
            let header = format!(">{node_id}\n");
            fasta_writer.write_all(header.as_bytes())?;
            let sequence = node.sequence.as_bytes();
            let sequence = crate::utils::insert_newlines(
                std::str::from_utf8(&sequence).unwrap(),
                120
            );
            fasta_writer.write_all(sequence.as_bytes())?;
            fasta_writer.write_all(b"\n")?;
        }
        
        Ok(())
    }

}

pub type GfaGraph = AsmGraph<String>;

