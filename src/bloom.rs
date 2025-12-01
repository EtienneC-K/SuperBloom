
//use bitvec::BitVec;
use seq_hash::NtHasher;
use bit_vec::BitVec;
use std::sync::Mutex;

//size of blocks, for now constants to fit a rather small L2 cache for labtops used in 2025
//pub const BLOCK_SIZE: usize = 1<<14; //2 097 152
//pub const NB_BLOCKS: usize = 1<<19; //16 384 for now, will see later to make it varaible

pub struct BloomFilter {
    pub filter: Vec<Mutex<BitVec>>,
    pub hashers: Vec<NtHasher>, //a vec of hash functions maybe ,or smth like an ntHash build je sais pas
    block_size: usize,
    nb_blocks: usize,
}

impl BloomFilter {
    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u32, k: usize, block_size: usize, 
        nb_blocks: usize) -> Self {

        let mut filter: Vec<Mutex<BitVec>> = Vec::new();
        for _ in 0..(size/block_size) {
            filter.push(Mutex::new(BitVec::from_elem(block_size, false)));
        }
        Self {
            //size,
            //n_hashes,
            //filter: vec![Mutex::new(BitVec::from_elem(BLOCK_SIZE, false)); size/BLOCK_SIZE],
            filter: filter,
            hashers: init_hashers(n_hashes, seed, k),
            block_size,
            nb_blocks,
        }
    }

    pub fn new(size: usize, n_hashes: usize, k: usize, block_size: usize, nb_blocks: usize) -> Self {
        let seed: u32 = 42;
        Self::new_with_seed(size, n_hashes, seed, k, block_size, nb_blocks)
    }

    ///checks if the kmer with specified minimizer hash, and multiple hashes is
    ///inside the bloom filter, inserts it if needed
    pub fn check_and_insert(&self, hashed_minimizer: u64, kmer_s_hashes: Vec<u32>) -> bool {
        let mut present: bool = true;
        let blocknum: usize = (hashed_minimizer as usize)%self.nb_blocks;
        let mut block = self.filter[blocknum].lock().unwrap();
        
        for hash in kmer_s_hashes {
            //to get the address, heavy bits are from the minimizer (giving the block)
            //and light bits are given by the hash of the kmer himself
            let address = hash as usize%self.block_size;
            if !block.get(address).unwrap_or(false) {
                block.set(address, true);
                present = false;
            }
        }
        present
    }

    //pub fn check_true_bits(&self) -> usize {
    //    let mut counter: usize = 0;
    //    for i in 0..4294967296 {
    //        if self.filter.get(i).unwrap() {
    //            counter += 1;
    //        }
    //    }
    //    counter
    //}
}


///to get the NtHasher hasher's when creating the bloomfilter
fn init_hashers(n_hashes : usize, seed: u32, k: usize) -> Vec<NtHasher> {
    let mut hasher_vec: Vec<NtHasher> = Vec::new();
    //we build hashers with slightly spaced seeds
    for i in 0..n_hashes {
        let hasher = <seq_hash::NtHasher>::new_with_seed(k, seed+(42*i as u32));
        hasher_vec.push(hasher);
    }
    hasher_vec
}
