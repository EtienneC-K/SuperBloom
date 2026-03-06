///module thatt utilizes some of imartayan's libraries to return minmizers as well as their
///positions, so they can later be hashed and their kmers un packed and hashed to be inserted

use crate::decyclers;

use packed_seq::{PackedSeqVec, SeqVec, Seq, PackedSeq};
use simd_minimizers::{canonical_minimizers};
use decyclers::{Decycler};

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

    (super_kmers, minimizer_vals, packed_seq) //for when ill use the rolling part of rolling
}

///takes same input and returns same thing as minimizers_x_positions plus decycling sets slice, 
///and uses those decycling sets to choose minimizers instead of random order
pub fn decycling_mins_x_pos (packed_seq: PackedSeqVec, k: u16, m: u16, decycler_set: &Decycler)
    -> (Vec<u32>, Vec<u64>, PackedSeqVec) {

    //protection against sequences too short
    if packed_seq.len() <= k as usize+2 {
        return ([].to_vec(), [].to_vec(), packed_seq);
    }

    let mut minimizer_vals = Vec::new();
    let mut super_kmers = Vec::new(); //those are actually positions, not the superkmers themselves

    //start by making an index to know if theyre in the decycling set
    let mut is_decycler: Vec<bool> = Vec::with_capacity(packed_seq.len()-m as usize+1);
    //let vec_ci = init_vec_ci(m);
    for i in 0..packed_seq.len()-m as usize+1 {
        is_decycler.push(decycler_set.lookup(packed_seq.slice(i..i+m as usize)));
    }

    let mut mini_addrs: Vec<usize> = Vec::with_capacity(packed_seq.len()-k as usize+1);
    //do a max and a is_decyc var, and look for the minimums in each chunk
    //first look for the minimizer in the first kmer
    let (mut min_addr, mut _is_decyc, mut _min_lexic) = 
        mins_from_kmer(packed_seq.as_slice(), &is_decycler, 0, m, k);
    mini_addrs.push(min_addr);

    //then we roll
    for i in 1..packed_seq.len()+1-k as usize {
        //we check if the old minimizer is still in the new window
        if min_addr < i {
            //is the previous memeber of decycling set rolls out, or it wasn't one in the first
            //place and it just rolled out, we look for a new minimizer
            (min_addr, _is_decyc, _min_lexic) =
                mins_from_kmer(packed_seq.as_slice(), &is_decycler, i, m, k);
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

    (super_kmers, minimizer_vals, packed_seq)
}

fn mins_from_kmer<'a>(packed_seq: PackedSeq<'a>, is_decycler: &Vec<bool>, i: usize, m: u16, k: u16) 
    -> (usize, bool, PackedSeq<'a>) {

    for j in i..i+k as usize-m as usize+1 {
        if is_decycler[j] {return (j, true, packed_seq.slice(j..j+m as usize))};
    }
    (i, false, packed_seq.slice(i..i+m as usize))
}
