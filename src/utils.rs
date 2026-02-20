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

///since no implementation directly on u128 exists, I just used whatever numbers of shifts
pub fn xorshift_u128(mut x: u128) -> u128 {
    x ^= x << 17;
    x ^= x >> 23;
    x ^= x << 5;
    x
}
