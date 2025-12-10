///module with small utility fucntions

//use bitvec::prelude::*;
//use packed_seq::{PackedSeq, Seq};

///hashing function for u64 using xorshift
pub fn xorshift_u64(mut x: u64) -> u64 {
    x ^= x<<13;
    x ^= x>>7;
    x ^= x<<17;
    x
}
