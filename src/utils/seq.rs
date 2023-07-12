use std::fmt;

#[derive(Debug)]
pub struct SeqInterval {
    pub tid: usize,
    pub beg: usize,
    pub end: usize,
}

impl fmt::Display for SeqInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.tid, self.beg, self.end)
    }
}