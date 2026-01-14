
//use bitvec::BitVec;
use seq_hash::NtHasher;
use bit_vec::BitVec;
use std::sync::Mutex;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::prelude::ParallelSliceMut;
use std::ops::Deref;

pub struct BloomFilter {
    pub filter: Vec<Mutex<Vec<BitVec>>>,
    pub hasher: NtHasher, //a vec of hash functions maybe ,or smth like an ntHash build je sais pas
    block_size: usize,
    nb_blocks: usize,
    pub n_hashes: usize,
}

impl BloomFilter {
    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u32, k: usize, block_size: usize, 
        nb_blocks: usize) -> Self {

        let magic_mutex_amount = 1024;
        assert!(magic_mutex_amount*block_size <= size);
        assert!(nb_blocks >= magic_mutex_amount);
        let mut filter: Vec<Mutex<Vec<BitVec>>> = Vec::new();
        for _ in 0..magic_mutex_amount {
            filter.push(Mutex::new(vec![BitVec::from_elem(block_size, false);
                                        size/(block_size*magic_mutex_amount)]));
        }
        Self {
            filter: filter,
            hasher: init_hasher(n_hashes, seed, k),
            block_size,
            nb_blocks,
            n_hashes,
        }
    }

    pub fn new(size: usize, n_hashes: usize, k: usize, block_size: usize, nb_blocks: usize) -> Self {
        let seed: u32 = 42;
        Self::new_with_seed(size, n_hashes, seed, k, block_size, nb_blocks)
    }

    ///checks if the kmer with specified minimizer hash, and multiple hashes is
    ///inside the bloom filter, inserts it if needed
    pub fn check_and_insert(&self, hashed_minimizer: u64, kmer_s_hashes: Vec<u64>) -> bool {
        let mut present: bool = true;
        let blocknum: usize = (hashed_minimizer as usize)%1024;
        let subblocknum: usize = ((hashed_minimizer as usize)/1024)%(self.nb_blocks/1024);
        let mut block = self.filter[blocknum].lock().unwrap();
        let mut subblock = &mut block[subblocknum];

        for hash in kmer_s_hashes {
            //to get the address, heavy bits are from the minimizer (giving the block)
            //and light bits are given by the hash of the kmer himself
            let address = hash as usize%self.block_size;
            if !subblock.get(address).unwrap() {
                subblock.set(address, true);
                present = false;
            }
        }
        present
    }

    ///now unusable and wrong because of change in format of the filter
    pub fn check_true_bits(&self) -> usize {
        //let mut counter: usize = 0;
        //for block in &self.filter {
        //    let unlocked_block = block.lock().unwrap();
        //    for i in 0..unlocked_block.len() {
        //        if unlocked_block.get(i).unwrap() {
        //            counter += 1;
        //        }
        //    }
        //}
        //counter
        1
    }

    ///counts different metrics like fill rate, avg fille rate of non empties, median one etc...
    ///returns : count of non epty blocks, max filled count, median filled count, avrg filled count
    pub fn count_it_all(&self) -> (usize, usize, usize, usize) {
        //first make a list with all non zero rates
        let counts_list: Mutex<Vec<usize>> = Mutex::new(Vec::new());
        let total_counter: Mutex<usize> = Mutex::new(0);
        let _ = &self.filter.iter().par_bridge().for_each(|block| {
            let unlocked_block = block.lock().unwrap(); //its a Vec<BitVec>
            for bit_vector in unlocked_block.deref() {
                let mut counter: usize = 0;
                for i in 0..bit_vector.len() {
                    if bit_vector.get(i).unwrap() {
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
            }
        });

        //once we have the list, its time to sort it
        let mut unlocked_counts_list = counts_list.lock().unwrap();
        unlocked_counts_list.par_sort_unstable();

        //now to calculate what we're looking for
        let non_zero_counters: usize = unlocked_counts_list.len();
        let max_counter: usize = unlocked_counts_list[unlocked_counts_list.len()-1];
        let median_counter: usize = unlocked_counts_list[unlocked_counts_list.len()/2 - 1];
        let average_counter: usize = *total_counter.lock().unwrap()/unlocked_counts_list.len();

        (non_zero_counters, max_counter, median_counter, average_counter)
    }
}


///to get the NtHasher hasher's when creating the bloomfilter
fn init_hasher(n_hashes : usize, seed: u32, k: usize) -> NtHasher {
    let mut hasher_vec: Vec<NtHasher> = Vec::new();
    //we build hashers with slightly spaced seeds
    let hasher = <seq_hash::NtHasher>::new_with_seed(k, seed);
    hasher
}
