///module to read an inputed fasta file, and maybe later a .txt file of file
///output will always be some compressed sequences (packed_seq) from imartayan

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use packed_seq::{PackedSeqVec, SeqVec};
use needletail::FastxReader;

//pub fn read_fasta(fasta_file: String) -> packed_seq::packed_seq::PackedSeqVecBase<2> {
pub fn _read_fasta(fasta_file: String) -> PackedSeqVec {
    //var that will contain the concatenation of all lines before conversion to packed_seq
    let mut full_ascii: String = String::new();
    if let Ok(lines) = _read_lines(fasta_file) {
        for line in lines {
            let unwrapped_line = line.expect("Problem reading a FASTA");
            let line_bytes = unwrapped_line.as_bytes();
            //filter out all the comments
            if line_bytes.len() >0 && line_bytes[0] != b'>' && line_bytes[0] != b';' {
                full_ascii += &unwrapped_line;
            }
        }
    }
    //once we have everything in a single String, its time to turn it into a more efficient
    //packed_seq, to be used quickly later
    let packed_seq = PackedSeqVec::from_ascii(full_ascii.as_bytes());

    packed_seq
}

///function for reading a file of file to handle lots of fasta at once
///the fof should have the path to a single fasta on each line
pub fn _read_fof(fof_file: String) -> Vec<String> {
    let mut iter_files: Vec<String> = Vec::new();
    if let Ok(lines) = _read_lines(fof_file) {
        for line in lines {
            iter_files.push(line.unwrap());
        }
    }
    iter_files
}

///classic function to simply read any file line by line efficiently
pub fn _read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

/// I never hated a struct more than an Iterator with lifetimes, worst invention in Mankind history
pub struct Hell {
    pub fxreader: Box<dyn FastxReader>,
    pub chunk_size: usize,
}

impl Iterator for Hell {
    type Item = Vec<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunk = Vec::new();
        for _ in 0..self.chunk_size {
            let result = match self.fxreader.next() {
                Some(res) => res,
                None => break,
            };
            let seq_red = result.unwrap().seq().to_mut().clone();
            chunk.push(seq_red);
        }
        if chunk.is_empty() {
            None
        } else {
            return Some(chunk);
        }
    }
}
