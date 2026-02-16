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

pub fn _xorshift_u32(mut x: u32) -> u32 {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    x
}

///computes the u64's that will be inserted in the block right after that
pub fn compute_insertions(relevant_addresses: &[usize]) -> Vec<u64> {
    let mut to_inserts: Vec<u64> = Vec::with_capacity(relevant_addresses.len());
    for address in relevant_addresses {
        to_inserts.push(1<<(63-address%64) as u64); //trust the calculation
    }
    to_inserts
}
