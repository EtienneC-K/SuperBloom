///module thatt utilizes some of imartayan's libraries to return minmizers as well as their
///positions, so they can later be hashed and their kmers un packed and hashed to be inserted

use crate::decyclers;

use std::cell::RefCell;
use packed_seq::{PackedSeqVec, SeqVec, Seq, PackedSeq};
use seq_hash::NtHasher;
use simd_minimizers::{canonical_minimizers};
use decyclers::{Decycler};

#[derive(Default)]
struct SimdMinimizerScratch {
    minimizer_positions: Vec<u32>,
}

thread_local! {
    static SIMD_MINIMIZER_SCRATCH: RefCell<SimdMinimizerScratch> =
        RefCell::new(SimdMinimizerScratch::default());
}

///function that does all the job we're looking for here, with given kmer and word lengths
///also converts the packedseqvec to a bitvec to easily slice it later on
pub fn minimizers_x_positions(packed_seq: PackedSeqVec, k: u16, m: u16) 
        -> (Vec<u32>, Vec<u64>, PackedSeqVec) {
    let max_windows = packed_seq.len().saturating_sub(k as usize) + 1;
    let mut super_kmers = Vec::with_capacity(max_windows); //those are actually positions, not the superkmers themselves
    let window_size = k-m+1;
    let minimizer_length = m;

    let minimizer_vals = SIMD_MINIMIZER_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        scratch.minimizer_positions.clear();
        let current_capacity = scratch.minimizer_positions.capacity();
        scratch
            .minimizer_positions
            .reserve(max_windows.saturating_sub(current_capacity));

        let seed: u32 = 42;
        let hasher = <NtHasher>::new_with_seed(minimizer_length.into(), seed);
        let output = canonical_minimizers(minimizer_length.into(), window_size.into())
            .hasher(&hasher)
            .super_kmers(&mut super_kmers)
            .run(packed_seq.as_slice(), &mut scratch.minimizer_positions);
        output.values_u64().collect()
    });

    (super_kmers, minimizer_vals, packed_seq)
}

///takes same input and returns same thing as minimizers_x_positions plus decycling sets slice, 
///and uses those decycling sets to choose minimizers instead of random order
pub fn decycling_mins_x_pos (packed_seq: PackedSeqVec, k: u16, m: u16, decycler_set: &Decycler)
    -> (Vec<u32>, Vec<u64>, PackedSeqVec) {

    //protection against sequences too short
    if packed_seq.len() < k as usize {
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

#[cfg(test)]
mod tests {
    use super::{decycling_mins_x_pos, minimizers_x_positions, mins_from_kmer};
    use crate::decyclers::Decycler;
    use packed_seq::{PackedSeqVec, Seq, SeqVec};

    #[test]
    fn simd_minimizers_return_aligned_positions_and_values() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTAC");
        let (super_kmers, minimizers, original) = minimizers_x_positions(sequence.clone(), 5, 3);

        assert_eq!(original.len(), sequence.len());
        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
        assert_eq!(super_kmers[0], 0);
        assert!(super_kmers.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn decycling_minimizers_return_empty_for_short_sequences() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACG");
        let (super_kmers, minimizers, original) = decycling_mins_x_pos(sequence.clone(), 4, 3, &decycler);

        assert!(super_kmers.is_empty());
        assert!(minimizers.is_empty());
        assert_eq!(original.len(), sequence.len());
    }

    #[test]
    fn decycling_minimizers_return_monotonic_super_kmer_starts() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTACGT");
        let (super_kmers, minimizers, _) = decycling_mins_x_pos(sequence, 5, 3, &decycler);

        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
        assert_eq!(super_kmers[0], 0);
        assert!(super_kmers.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn mins_from_kmer_prefers_first_decycling_candidate() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTAC");
        let is_decycler = vec![false, true, true, false];
        let (addr, is_member, minimizer) = mins_from_kmer(sequence.as_slice(), &is_decycler, 0, 3, 5);

        assert_eq!(addr, 1);
        assert!(is_member);
        assert_eq!(minimizer.as_u64(), sequence.slice(1..4).as_u64());
    }

    #[test]
    fn mins_from_kmer_falls_back_to_window_start_when_needed() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTAC");
        let is_decycler = vec![false, false, false, false];
        let (addr, is_member, minimizer) = mins_from_kmer(sequence.as_slice(), &is_decycler, 0, 3, 5);

        assert_eq!(addr, 0);
        assert!(!is_member);
        assert_eq!(minimizer.as_u64(), sequence.slice(0..3).as_u64());
    }

    #[test]
    fn mins_from_kmer_ignores_candidates_outside_window() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTAC");
        let is_decycler = vec![false, false, false, true];
        let (addr, is_member, minimizer) = mins_from_kmer(sequence.as_slice(), &is_decycler, 0, 3, 5);

        assert_eq!(addr, 0);
        assert!(!is_member);
        assert_eq!(minimizer.as_u64(), sequence.slice(0..3).as_u64());
    }

    #[test]
    fn simd_minimizers_support_exact_kmer_length() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTA");
        let (super_kmers, minimizers, _) = minimizers_x_positions(sequence, 5, 3);

        assert_eq!(super_kmers.len(), 1);
        assert_eq!(minimizers.len(), 1);
        assert_eq!(super_kmers[0], 0);
    }

    #[test]
    fn simd_minimizers_return_original_sequence_unchanged() {
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let expected = sequence.as_slice().as_u64();
        let (_, _, original) = minimizers_x_positions(sequence, 5, 3);

        assert_eq!(original.as_slice().as_u64(), expected);
    }

    #[test]
    fn decycling_minimizers_become_non_empty_at_k_plus_three() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACG");
        let (super_kmers, minimizers, _) = decycling_mins_x_pos(sequence, 4, 3, &decycler);

        assert!(!super_kmers.is_empty());
        assert_eq!(super_kmers.len(), minimizers.len());
    }

    #[test]
    fn decycling_minimizer_values_fit_requested_mmer_width() {
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGTACGT");
        let (super_kmers, minimizers, original) = decycling_mins_x_pos(sequence, 5, 3, &decycler);

        assert_eq!(super_kmers.len(), minimizers.len());
        for minimizer in minimizers {
            assert!(minimizer < (1 << 6));
        }
        assert_eq!(original.len(), 12);
    }
}
