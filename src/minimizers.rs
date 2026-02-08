///module thatt utilizes some of imartayan's libraries to return minmizers as well as their
///positions, so they can later be hashed and their kmers un packed and hashed to be inserted

use crate::decyclers;

use packed_seq::{PackedSeqVec, SeqVec, Seq};
use simd_minimizers::{canonical_minimizers};
use decyclers::Decycler;
//use bitvec::prelude::*;
//use seq_hash::{NtHasher};

///function that does all the job we're looking for here, with given kmer and word lengths
///also converts the packedseqvec to a bitvec to easily slice it later on
pub fn minimizers_x_positions(packed_seq: PackedSeqVec, k: u16, m: u16) 
        -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    //variables that need to be initialized for the canonical minimizer to be computed
    let mut minimizer_positions = Vec::new();
    let mut super_kmers = Vec::new(); //those are actually positions, not the superkmers themselves
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

///takes same input and returns same thing as minimizers_x_positions plus decycling sets slice, 
///and uses those decycling sets to choose minimizers instead of random order
pub fn decycling_mins_x_pos (packed_seq: PackedSeqVec, k: u16, m: u16, decycler_set: &Decycler)
    -> (Vec<u32>, Vec<u64>, PackedSeqVec) {

    let mut minimizer_vals = Vec::new();
    let mut super_kmers = Vec::new(); //those are actually positions, not the superkmers themselves

    //start by making an index to know if theyre in the decycling set
    let mut is_decycler: Vec<bool> = Vec::with_capacity(packed_seq.len()-m as usize+1);
    for i in 0..packed_seq.len()-m as usize+1 {
        is_decycler.push(decycler_set.lookup(packed_seq.slice(i..i+m as usize)));
    }

    let mut mini_addrs: Vec<usize> = Vec::with_capacity(packed_seq.len()-k as usize+1);
    //do a max and a is_decyc var, and look for the minimums in each chunk
    for i in 0..packed_seq.len()+1-k as usize {
        //find the minimizer in this kmer
        let mut is_decyc: bool = is_decycler[i];
        let mut min_lexic = packed_seq.slice(i..i+m as usize);
        let mut min_addr: usize = i;
        for j in i+1..i+k as usize-m as usize+1 {
            let current_decyc: bool = is_decycler[j];
            let current_lexic = packed_seq.slice(j..j+m as usize);
            if current_decyc && !is_decyc {
                //we fond a first member of a decycling set
                is_decyc = true;
                min_lexic = current_lexic;
                min_addr = j;
            } else if is_decyc == current_decyc && min_lexic > current_lexic {
                //found a lexicographically smaller one
                min_lexic = current_lexic;
                min_addr = j;
            }
        }
        mini_addrs.push(min_addr);
    }

    //using the list of the minimizer positions, build the super_kmers
    let mut first_superkmer_addr: u32 = 0;
    let mut last_minimizer_pos: usize = mini_addrs[0];
    for i in 1..packed_seq.len()+1-k as usize {
        if mini_addrs[i] != last_minimizer_pos {
            //on passe au superkmer suivant
            last_minimizer_pos = mini_addrs[i];
            //let last_superkmer_addr: usize = i-1;
            super_kmers.push(first_superkmer_addr);
            minimizer_vals.push(
                packed_seq.slice(mini_addrs[first_superkmer_addr as usize]
                ..mini_addrs[first_superkmer_addr as usize]+m as usize)
                .as_u64());
            first_superkmer_addr = i as u32;
        }
    }
    //not forgetting to close the last superkmer
    super_kmers.push(first_superkmer_addr);
    minimizer_vals.push(
        packed_seq.slice(mini_addrs[first_superkmer_addr as usize]
        ..mini_addrs[first_superkmer_addr as usize]+m as usize)
        .as_u64());

    //once we have the startign positions of all the super kmers
    //we can get the address, and then value through slices as u64 of their minimizers
    //for super_kmer in super_kmers {
    //    minimizer_vals.push(
    //        packed_seq.slice(mini_addrs[super_kmer as usize]
    //        ..mini_addrs[super_kmer as usize]+m as usize)
    //        .as_u64());
    //}

    (super_kmers, minimizer_vals, packed_seq)
}
