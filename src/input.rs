///module to read an inputed fasta file, and maybe later a .txt file of file
///output will always be some compressed sequences (packed_seq) from imartayan

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
//use std::str;
use packed_seq::{PackedSeqVec, SeqVec};
//use bitvec::prelude::*;
use async_stream::stream;
use futures_core::stream::Stream;

//pub fn read_fasta(fasta_file: String) -> packed_seq::packed_seq::PackedSeqVecBase<2> {
pub fn read_fasta(fasta_file: String) -> PackedSeqVec {
    //var that will contain the concatenation of all lines before conversion to packed_seq
    let mut full_ascii: String = String::new();
    if let Ok(lines) = read_lines(fasta_file) {
        for line in lines {
            let unwrapped_line = line.expect("Problem reading a FASTA");
            let line_bytes = unwrapped_line.as_bytes();
            //filter out all the comments
            if line_bytes[0] != b'>' && line_bytes[0] != b';' {
                full_ascii += &unwrapped_line;
            }
        }
    }
    //once we have everything in a single String, its time to turn it into a more efficient
    //packed_seq, to be used quickly later
    let packed_seq = PackedSeqVec::from_ascii(full_ascii.as_bytes());

    //BitVec that i made, uses same encoding as imartayan but with a type that i know better
    /*let mut bitvec_seq = bitvec![0; full_ascii.len()*2];
    for (i, letter) in full_ascii.enumerate() {
        if letter == "C" {
            bitvec_seq[2*i+1] = 1;
        } else if letter == "T" {
            bitvec_seq[2*i] = 1;
        } else if letter == "G" {
            bitvec_seq[2*i+1] = 1;
            bitvec_seq[2*i] = 1;
        }
    }
    packed_seq, bitvec_seq*/
    packed_seq
}

///function for reading a file of file to handle lots of fasta at once
///the fof should have the path to a single fasta on each line
pub fn read_fof(fof_file: String) -> impl Stream<Item = PackedSeqVec> {
    if let Ok(lines) = read_lines(fof_file) {
        stream! {
            for line in lines {
                let unwrapped_line = line.expect("Problem reading fof");
                yield(read_fasta(unwrapped_line));
            }
        }
    } else {
        //sinon on renvoie le truc le plus lame possible
        panic!("Des bigs problèmes sur la fonction de lecture de fof");
    }
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
