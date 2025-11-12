///module with small utility fucntions

use bitvec::prelude::*;
use packed_seq::{PackedSeq, Seq};

///hashing function for u64 using xorshift
pub fn xorshift_u64(mut x: u64) -> u64 {
    x ^= x<<13;
    x ^= x>>7;
    x ^= x<<17;
    x
}

///converts the packedseq sequence to a BitVec for manipulation in the hash table
///table and such
pub fn convert_seqkmer(sequence: PackedSeq) -> BitVec {
    //sequence.as_slice().iter_bp().enumerate() //hopefully it works without slice
    //let mut bit_seq = bitvec![0; 2*sequence.len()];
    let mut bit_seq = bitvec![0; 64]; //for now fixed size 64 for encoding of any kmer
    for (i, nuc) in sequence.iter_bp().enumerate() {
        if nuc == 1 {
            bit_seq.set(2*i+1, true);
        } else if nuc == 2 {
            bit_seq.set(2*i, true);
        } else if nuc == 3 {
            bit_seq.set(2*i+1, true);
            bit_seq.set(2*i, true);
        }
    }
    bit_seq
}
