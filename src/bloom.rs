use crate::super_bitvec;
use crate::minimizers;
use crate::decyclers;
use crate::utils;

use seq_hash::NtHasher;
use bit_vec::BitVec;
use super_bitvec::SuperBitVec;
use std::sync::Mutex;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::prelude::ParallelSliceMut;
use std::ops::Deref;
use packed_seq::{PackedSeqVec, SeqVec, PackedSeq, Seq};
use decyclers::Decycler;
use utils::{xorshift_u64, xorshift_u128, sum_vec_bool};
use minimizers::{decycling_mins_x_pos, minimizers_x_positions};
use rand::prelude::*;

pub struct BloomFilter {
    pub filter: Vec<Mutex<Vec<SuperBitVec>>>,
    pub block_size: usize,
    pub nb_blocks: usize,
    pub n_hashes: usize,
    block_size_mask: usize,
}

impl BloomFilter {
    pub fn new_with_seed(size: usize, n_hashes: usize, _seed: u32, _k: usize, block_size: usize, 
        nb_blocks: usize) -> Self {

        let magic_mutex_amount = 1024;
        assert!(magic_mutex_amount*block_size <= size);
        assert!(nb_blocks >= magic_mutex_amount);
        let mut filter: Vec<Mutex<Vec<SuperBitVec>>> = Vec::new();
        for _ in 0..magic_mutex_amount {
            filter.push(Mutex::new(vec![SuperBitVec::new(block_size);
                                        size/(block_size*magic_mutex_amount)]));
        }
        Self {
            filter: filter,
            block_size,
            nb_blocks,
            n_hashes,
            block_size_mask: block_size-1,
        }
    }

    pub fn new(size: usize, n_hashes: usize, k: usize, block_size: usize, nb_blocks: usize) -> Self {
        let seed: u32 = 42;
        Self::new_with_seed(size, n_hashes, seed, k, block_size, nb_blocks)
    }

    ///checks if the kmer with specified minimizer hash, and multiple hashes is
    ///inside the bloom filter, inserts it if needed
    pub fn _check_and_insert(&self, subblock: &mut BitVec, mut hash: u64) -> bool {
        let mut present: bool = true;

        for _i in 0..self.n_hashes {
            let address = hash as usize%self.block_size;
            if !subblock.get(address).unwrap() {
                subblock.set(address, true);
                present = false;
            }
            hash = xorshift_u64(hash);
        }
        present
    }


    ///counts different metrics like fill rate, avg fille rate of non empties, median one etc...
    ///returns : count of non empty blocks, max filled count, median filled count, avrg filled count
    pub fn count_it_all(&self) -> (usize, usize, usize, usize, usize) {
        //first make a list with all non zero rates
        let counts_list: Mutex<Vec<usize>> = Mutex::new(Vec::new());
        let total_counter: Mutex<usize> = Mutex::new(0);
        let filled_counter: Mutex<usize> = Mutex::new(0);
        let _ = &self.filter.iter().par_bridge().for_each(|block| {
            let unlocked_block = block.lock().unwrap(); //its a Vec<BitVec>
            for bit_vector in unlocked_block.deref() {
                let mut counter: usize = 0;
                for i in 0..bit_vector.len() {
                    if bit_vector.get(i) {
                        counter += 1;
                    }
                }
                if counter > 0 {
                    let mut el_liste = counts_list.lock().unwrap();
                    el_liste.push(counter);
                    drop(el_liste);
                    let mut el_counter = total_counter.lock().unwrap();
                    *el_counter = el_counter.saturating_add(counter);
                    drop(el_counter)
                }

                //add a counter of blocks filled above a certain threshhold
                let threshhold: f64 = 0.9;
                if counter as f64/self.block_size as f64 > threshhold {
                    let mut el_filled_counter = filled_counter.lock().unwrap();
                    *el_filled_counter = el_filled_counter.saturating_add(1);
                    drop(el_filled_counter);
                }
            }
        });

        //once we have the list, its time to sort it
        let mut unlocked_counts_list = counts_list.lock().unwrap();
        if unlocked_counts_list.is_empty() {
            return (0, 0, 0, 0, 0);
        }
        unlocked_counts_list.par_sort_unstable();

        //now to calculate what we're looking for
        let non_zero_counters: usize = unlocked_counts_list.len();
        let max_counter: usize = unlocked_counts_list[unlocked_counts_list.len()-1];
        let median_counter: usize = unlocked_counts_list[unlocked_counts_list.len()/2 - 1];
        let average_counter: usize = *total_counter.lock().unwrap()/unlocked_counts_list.len();
        let filled_count = *filled_counter.lock().unwrap();

        (non_zero_counters, max_counter, median_counter, average_counter, filled_count)
    }

    ///checks the false negative and false positive counts of the bloom filter
    pub fn count_false_bloom(&self, to_check: Vec<PackedSeqVec>, k: u16, m: u16, l: u16, decycler_set: &Decycler) -> (f64, f64) {
        let (false_negs, total_neg_tests, nb_seq_neg_tests) = self.count_false_negatives(to_check, k, m, l, decycler_set);
        let false_pos = self.count_false_positives(k, m, l, total_neg_tests, nb_seq_neg_tests, decycler_set);
        (false_negs, false_pos)
    }


    ///using a set of kmer that where supposed to be inserted, and randomly generated kmers checks
    ///the rates of false negatives 
    pub fn count_false_negatives (
        &self,
        to_check : Vec<PackedSeqVec>,
        k: u16,
        m: u16,
        l: u16,
        decycler_set: &Decycler,
        ) -> (f64, usize, usize) {
        //start by checking for false negatives
        let nb_seq_neg_tests: usize = to_check.len();
        let mut false_negative_count: usize = 0;
        let mut total_count: usize = 0;
        for sequence in to_check {
            total_count += sequence.len()-(k as usize)+1;
            let (_count_true, count_false): (usize, usize);
            if l <= 31 {
                let presence_vec = self.check_sequence(sequence, k, m, l, decycler_set);
                count_false = presence_vec.len()-sum_vec_bool(&presence_vec);
            } else {
                let presence_vec = self.check_sequence_u128(sequence, k, m, l, decycler_set);
                count_false = presence_vec.len()-sum_vec_bool(&presence_vec);
            }
            false_negative_count += count_false;
        }
        let false_proportion: f64 = false_negative_count as f64/total_count as f64;
        (false_proportion, total_count, nb_seq_neg_tests)
    }

    ///generates random kmers, that are therefore likely not supposed to be here, and counts how
    ///many return as positive
    pub fn count_false_positives(&self, k: u16, m: u16, l: u16,
        total_false_negs: usize,
        nb_sequence_false_negs: usize,
        decycler_set: &Decycler,
        ) -> f64 {

        //protection against small examples with no good dice rolls
        if nb_sequence_false_negs < 1 {
            return -1.0;
        }

        let mut total_false_pos: usize = 0;
        let avg_len: usize = total_false_negs/nb_sequence_false_negs;
        for _i in 0..nb_sequence_false_negs {
            total_false_pos += self.make_n_check_sequence(k, m, l, avg_len, decycler_set);
        }

        let false_pos_rate: f64 = 
            total_false_pos as f64/((avg_len-k as usize)*nb_sequence_false_negs) as f64;
        false_pos_rate
    }

    ///generates a random sequence before counting false positives in it
    fn make_n_check_sequence(&self, k: u16, m: u16, l: u16, avg_len: usize, decycler_set: &Decycler) -> usize {
        //generate random sequence
        let mut rng = rand::rng();
        let mut seq: String = String::from("");
        let rand_mapping: String = String::from("ACGT");
        for _i in 0..avg_len {
            let mut rand_num: u8 = rng.random::<u8>();
            rand_num %= 4;
            seq.push(rand_mapping.as_bytes()[rand_num as usize] as char);
        }
        let sequence: PackedSeqVec = PackedSeqVec::from_ascii(seq.as_bytes());

        //now to check false positives
        let (count_true, _count_false): (usize, usize);
        if l <= 31 {
            let presence_vec = self.check_sequence(sequence, k, m, l, decycler_set);
            count_true = sum_vec_bool(&presence_vec);
        } else {
            let presence_vec = self.check_sequence_u128(sequence, k, m, l, decycler_set);
            count_true = sum_vec_bool(&presence_vec);
        }

        count_true
    }

    ///simply checks if a sequence of kmer is present or not, does no insertion and isn't thought to be suited
    ///for parallel operations, as only small checks at the end
    ///returns a vec of boolean with the i-th indexed boolean corresponding to if the i-th kmer
    ///returns a positive
    pub fn check_sequence(&self, original_sequence: PackedSeqVec, k: u16, m: u16, s: u16, decycler_set: &Decycler) -> Vec<bool> {
        //protection against small sequences
        if original_sequence.len() < k as usize {
            return vec![];
        }

        let mut presence_vec: Vec<bool> = Vec::with_capacity(original_sequence.len()-k as usize+1);
        let address_mask = (self.nb_blocks-1)>>10;


        let (super_kmers_positions, minimizers, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec);
        if decycler_set.m > 1 {
            //we use the decycler
            (super_kmers_positions, minimizers, sequence) =
                decycling_mins_x_pos(original_sequence, k, m, decycler_set);
        } else {
            //we use simd_minimizer
            (super_kmers_positions, minimizers, sequence) =
                minimizers_x_positions(original_sequence, k, m);
        }


        for i in 0..super_kmers_positions.len() {
            let hashed_minimizer: u64 = xorshift_u64(minimizers[i]);
            let start_pos: usize = super_kmers_positions[i] as usize;
            let end_pos: usize = if i==super_kmers_positions.len()-1 {sequence.len()+1-k as usize} 
                                    else {super_kmers_positions[i+1] as usize};
            //must compute the subblock by ourselves, its not furnished this time around
            let blocknum: usize = (hashed_minimizer as usize)%1024;
            let subblocknum: usize = ((hashed_minimizer as usize)>>10)&address_mask;
            let mut block = self.filter[blocknum].lock().unwrap();
            let subblock = &mut block[subblocknum];

            for j in start_pos..end_pos {
                let kmer: PackedSeq = sequence.slice(j..j+k as usize);
                let present: bool = self.check_kmer(subblock, kmer, s);
                presence_vec.push(present);
            }
            drop(block);
        }
        presence_vec
    }

    ///checks if a kmer is present
    fn check_kmer(&self, subblock: &mut SuperBitVec, kmer: PackedSeq, s: u16) -> bool {

        for i in 0..kmer.len()-s as usize +1 {
            let smer = kmer.slice(i..i+s as usize);
            let mut hash = xorshift_u64(smer.as_u64());
            for _j in 0..self.n_hashes {
                let address = hash as usize&self.block_size_mask;
                if !subblock.get(address) {
                    return false
                }
                hash = xorshift_u64(hash);
            }
        }
        true
    }

    pub fn check_sequence_u128(&self, original_sequence: PackedSeqVec, k: u16, m: u16, s: u16, decycler_set: &Decycler) -> Vec<bool> {
        //protection against small sequences
        if original_sequence.len() < k as usize {
            return vec![];
        }

        let mut presence_vec: Vec<bool> = Vec::with_capacity(original_sequence.len()-k as usize+1);
        let address_mask = (self.nb_blocks-1)>>10;


        let (super_kmers_positions, minimizers, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec);
        if decycler_set.m > 1 {
            //we use the decycler
            (super_kmers_positions, minimizers, sequence) =
                decycling_mins_x_pos(original_sequence, k, m, decycler_set);
        } else {
            //we use simd_minimizer
            (super_kmers_positions, minimizers, sequence) =
                minimizers_x_positions(original_sequence, k, m);
        }


        for i in 0..super_kmers_positions.len() {
            let hashed_minimizer: u64 = xorshift_u64(minimizers[i]);
            let start_pos: usize = super_kmers_positions[i] as usize;
            let end_pos: usize = if i==super_kmers_positions.len()-1 {sequence.len()+1-k as usize} 
                                    else {super_kmers_positions[i+1] as usize};
            //must compute the subblock by ourselves, its not furnished this time around
            let blocknum: usize = (hashed_minimizer as usize)%1024;
            let subblocknum: usize = ((hashed_minimizer as usize)>>10)&address_mask;
            let mut block = self.filter[blocknum].lock().unwrap();
            let subblock = &mut block[subblocknum];

            for j in start_pos..end_pos {
                let kmer: PackedSeq = sequence.slice(j..j+k as usize);
                let present: bool = self.check_kmer_u128(subblock, kmer, s);
                presence_vec.push(present);
            }
            drop(block);
        }
        presence_vec
    }

    ///checks if a kmer is present
    fn check_kmer_u128(&self, subblock: &mut SuperBitVec, kmer: PackedSeq, s: u16) -> bool {

        for i in 0..kmer.len()-s as usize+1 {
            let smer = kmer.slice(i..i+s as usize);
            let mut hash = xorshift_u128(smer.as_u128());
            for _j in 0..self.n_hashes {
                let address = hash as usize&self.block_size_mask;
                if !subblock.get(address) {
                    return false
                }
                hash = xorshift_u128(hash);
            }
        }
        true
    }



}


///to get the NtHasher hasher's when creating the bloomfilter
fn _init_hasher(seed: u32, k: usize) -> NtHasher {
    //we build hashers with different seeds
    let hasher = <seq_hash::NtHasher>::new_with_seed(k, seed);
    hasher
}

#[cfg(test)]
mod tests {
    use super::BloomFilter;
    use crate::decyclers::Decycler;
    use crate::super_bitvec::SuperBitVec;
    use crate::utils::xorshift_u64;
    use packed_seq::{PackedSeqVec, Seq, SeqVec};

    fn insert_sequence_with_decycler(
        bloom: &BloomFilter,
        sequence: &PackedSeqVec,
        k: u16,
        m: u16,
        s: u16,
        decycler: &Decycler,
    ) {
        let address_mask = bloom.block_size - 1;
        let nb_blocks = bloom.nb_blocks;
        let (super_kmers_positions, minimizer_values, packed) =
            crate::minimizers::decycling_mins_x_pos(sequence.clone(), k, m, decycler);

        for i in 0..super_kmers_positions.len() {
            let start = super_kmers_positions[i] as usize;
            let end = if i == super_kmers_positions.len() - 1 {
                packed.len() + 1 - k as usize
            } else {
                super_kmers_positions[i + 1] as usize
            };
            let hashed_minimizer = xorshift_u64(minimizer_values[i]) & (nb_blocks as u64 - 1);
            let blocknum = (hashed_minimizer as usize) & 1023;
            let subblocknum = ((hashed_minimizer as usize) >> 10) & ((bloom.nb_blocks >> 10) - 1);
            let mut block = bloom.filter[blocknum].lock().unwrap();
            let subblock = &mut block[subblocknum];

            for j in start..end {
                let kmer = packed.slice(j..j + k as usize);
                for offset in 0..kmer.len() - s as usize + 1 {
                    let smer = kmer.slice(offset..offset + s as usize);
                    if s <= 31 {
                        let mut hash = xorshift_u64(smer.as_u64());
                        for _ in 0..bloom.n_hashes {
                            let address = hash as usize & address_mask;
                            if !subblock.get(address) {
                                subblock.set(address, true);
                            }
                            hash = xorshift_u64(hash);
                        }
                    } else {
                        let mut hash = crate::utils::xorshift_u128(smer.as_u128());
                        for _ in 0..bloom.n_hashes {
                            let address = hash as usize & address_mask;
                            if !subblock.get(address) {
                                subblock.set(address, true);
                            }
                            hash = crate::utils::xorshift_u128(hash);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn check_and_insert_is_idempotent() {
        let bloom = BloomFilter::new(1 << 20, 3, 31, 1 << 10, 1 << 10);
        let mut subblock = bit_vec::BitVec::from_elem(1 << 10, false);
        let hash = 123456789_u64;

        assert!(!bloom._check_and_insert(&mut subblock, hash));
        assert!(bloom._check_and_insert(&mut subblock, hash));
    }

    #[test]
    fn check_sequence_returns_empty_for_short_sequences() {
        let bloom = BloomFilter::new(1 << 20, 3, 31, 1 << 10, 1 << 10);
        let decycler = Decycler::new(1);
        let sequence = PackedSeqVec::from_ascii(b"ACG");

        assert!(bloom.check_sequence(sequence, 5, 3, 3, &decycler).is_empty());
    }

    #[test]
    fn inserted_sequence_is_reported_present() {
        let bloom = BloomFilter::new(1 << 20, 3, 5, 1 << 10, 1 << 10);
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");

        insert_sequence_with_decycler(&bloom, &sequence, 5, 3, 3, &decycler);
        let results = bloom.check_sequence(sequence, 5, 3, 3, &decycler);

        assert_eq!(results.len(), 4);
        assert!(results.iter().all(|present| *present));
    }

    #[test]
    fn inserted_sequence_is_reported_present_for_u128_path() {
        let bloom = BloomFilter::new(1 << 20, 3, 40, 1 << 10, 1 << 10);
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );

        insert_sequence_with_decycler(&bloom, &sequence, 40, 3, 32, &decycler);
        let results = bloom.check_sequence_u128(sequence, 40, 3, 32, &decycler);

        assert_eq!(results.len(), 9);
        assert!(results.iter().all(|present| *present));
    }

    #[test]
    #[should_panic]
    fn bloom_new_panics_when_size_is_too_small() {
        let _ = BloomFilter::new(1024, 3, 31, 2, 1024);
    }

    #[test]
    #[should_panic]
    fn bloom_new_panics_when_nb_blocks_is_too_small() {
        let _ = BloomFilter::new(1 << 20, 3, 31, 1 << 10, 512);
    }

    #[test]
    fn count_it_all_is_zero_for_empty_bloom() {
        let bloom = BloomFilter::new(1 << 20, 3, 31, 1 << 10, 1 << 10);
        assert_eq!(bloom.count_it_all(), (0, 0, 0, 0, 0));
    }

    #[test]
    fn count_it_all_reports_non_zero_after_insertion() {
        let bloom = BloomFilter::new(1 << 20, 3, 5, 1 << 10, 1 << 10);
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");

        insert_sequence_with_decycler(&bloom, &sequence, 5, 3, 3, &decycler);
        let (non_zero, max_fill, median_fill, average_fill, _) = bloom.count_it_all();

        assert!(non_zero > 0);
        assert!(max_fill > 0);
        assert!(median_fill > 0);
        assert!(average_fill > 0);
    }

    #[test]
    fn check_kmer_returns_false_on_empty_subblock() {
        let bloom = BloomFilter::new(1 << 20, 3, 5, 1 << 10, 1 << 10);
        let mut subblock = SuperBitVec::new(1 << 10);
        let kmer = PackedSeqVec::from_ascii(b"ACGTA");

        assert!(!bloom.check_kmer(&mut subblock, kmer.as_slice(), 3));
    }

    #[test]
    fn check_kmer_u128_returns_false_on_empty_subblock() {
        let bloom = BloomFilter::new(1 << 20, 3, 40, 1 << 10, 1 << 10);
        let mut subblock = SuperBitVec::new(1 << 10);
        let kmer = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );

        assert!(!bloom.check_kmer_u128(&mut subblock, kmer.as_slice(), 32));
    }

    #[test]
    fn exact_kmer_sequence_is_reported_present() {
        let bloom = BloomFilter::new(1 << 20, 3, 5, 1 << 10, 1 << 10);
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(b"ACGTA");

        insert_sequence_with_decycler(&bloom, &sequence, 5, 3, 3, &decycler);
        let results = bloom.check_sequence(sequence, 5, 3, 3, &decycler);

        assert_eq!(results, vec![true]);
    }

    #[test]
    fn exact_kmer_sequence_is_reported_present_for_u128_path() {
        let bloom = BloomFilter::new(1 << 20, 3, 40, 1 << 10, 1 << 10);
        let mut decycler = Decycler::new(3);
        decycler.compute_blocks();
        let sequence = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );

        insert_sequence_with_decycler(&bloom, &sequence, 40, 3, 32, &decycler);
        let results = bloom.check_sequence_u128(sequence, 40, 3, 32, &decycler);

        assert_eq!(results, vec![true]);
    }

    #[test]
    fn count_false_positives_returns_minus_one_without_negative_samples() {
        let bloom = BloomFilter::new(1 << 20, 3, 31, 1 << 10, 1 << 10);
        let decycler = Decycler::new(1);

        assert_eq!(bloom.count_false_positives(5, 3, 3, 0, 0, &decycler), -1.0);
    }
}
