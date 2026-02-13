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

///makes pre_insertions in a local vecotr before locking the mutex and doing all the OR's
pub fn pre_insertions(thread_local_vec: &mut Vec<u64>, relevant_addresses: &[usize]) {
    //start by resetting the local vector to 0
    thread_local_vec.fill(0);

    //do all the small little or's like before
    for address in relevant_addresses {
        local_insert(thread_local_vec, *address);
    }
}

//TODO
//performs the creation of a u64 and inserts it in the local vector
fn local_insert(thread_local_vec: &mut Vec<u64>, addr: usize) {
    let to_insert: u64 = 1<<(63-addr%64) as u64;
    let block_num: usize = addr/64;
    thread_local_vec[block_num] = thread_local_vec[block_num] | to_insert;
}
