///module thatt utilizes some of imartayan's libraries to return minmizers as well as their
///positions, so they can later be hashed and their kmers un packed and hashed to be inserted

use packed_seq::{PackedSeqVec, SeqVec};
use simd_minimizers::{canonical_minimizers};
//use bitvec::prelude::*;
//use seq_hash::{NtHasher};

///function that does all the job we're looking for here, with given kmer and word lengths
///also converts the packedseqvec to a bitvec to easily slice it later on
pub fn minimizers_x_positions(packed_seq: PackedSeqVec, k: u16, m: u16) 
        -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    //variables that need to be initialized for the canonical minimizer to be computed
    let mut minimizer_positions = Vec::new();
    let mut super_kmers = Vec::new();
    //ce que je veux c'est super_kmers et les canonical_minimizers en vrai
    let window_size = k-m+1;
    let minimizer_length = m;

    //hasher
    let seed: u32 = 42;
    let hasher = <seq_hash::NtHasher>::new_with_seed(minimizer_length.into(), seed); //static seed 42 for now
    
    //actual computation from the library
    let minimizer_vals: Vec<u64> = canonical_minimizers(minimizer_length.into(), window_size.into())
        .hasher(&hasher)
        .super_kmers(&mut super_kmers)
        .run(packed_seq.as_slice(), &mut minimizer_positions)
        .values_u64()
        .collect();

    //(super_kmers, minimizer_vals)
    (super_kmers, minimizer_vals, packed_seq) //for when ill use the rolling part of rolling
    //hashes again, but not yet
}
