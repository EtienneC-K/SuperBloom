///module to read an inputed fasta file, and maybe later a .txt file of file
///output will always be some compressed sequences (packed_seq) from imartayan

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
//use std::str;
use packed_seq::{PackedSeqVec, SeqVec};
//use bitvec::prelude::*;
use needletail::FastxReader;
use needletail::parser::SequenceRecord;
use needletail::errors::ParseError;

//pub fn read_fasta(fasta_file: String) -> packed_seq::packed_seq::PackedSeqVecBase<2> {
pub fn read_fasta(fasta_file: String) -> PackedSeqVec {
    //var that will contain the concatenation of all lines before conversion to packed_seq
    let mut full_ascii: String = String::new();
    if let Ok(lines) = read_lines(fasta_file) {
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
pub fn read_fof(fof_file: String) -> Vec<String> {
    let mut iter_files: Vec<String> = Vec::new();
    if let Ok(lines) = read_lines(fof_file) {
        for line in lines {
            iter_files.push(line.unwrap());
        }
    }
    iter_files
}

//function that reads a fasta file with multiple (really many) reads, usefull for having plenty of
//reads without paying the header cost of having millions of files
//pub fn read_multi_fasta(fasta_file: String) -> 
//     maybe later idk                               //

///classic function to simply read any file line by line efficiently
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

//pub struct FastaReader<T: FastxReader> {
//    pub lines: T,
//    pub chunk_size: usize,
//}
//
//impl<'a, T> Iterator for FastaReader<T> {
//        type Item<'b> = Vec<Option<Result<SequenceRecord<'b>, ParseError>>>;
//
//        fn next(&mut self) -> Option<Self::Item> {
//            let mut chunk = Vec::with_capacity(self.chunk_size);
//            for _ in 0..self.chunk_size {
//                if let Some(line) = self.lines.next() {
//                    match line {
//                        Ok(content) => chunk.push(line),
//                        Err(_) => break,
//                    }
//                } else {
//                    break;
//                }
//            }
//            if chunk.is_empty() {
//                None
//            } else {
//                Some(chunk)
//            }
//        }
//}

/// I never hated a struct more than an Iterator with lifetimes, worst invention in Mankind history
pub struct Hell {
    pub fxreader: Box<dyn FastxReader>,
    pub chunk_size: usize,
}

impl Iterator for Hell {
    type Item = Vec<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.fxreader.next()?;
        let seq_red = result.unwrap().seq().to_mut().clone();
        let singleton_chunk = vec![seq_red];
        return Some(singleton_chunk);
    }
}
