
//use bitvec::BitVec;
use seq_hash::NtHasher;
use bit_vec::BitVec;

//size of blocks, for now constants to fit a rather small L2 cache for labtops used in 2025
pub const BLOCK_SIZE: usize = 1<<21; //2 097 152

pub struct BloomFilter {
    //size: usize,
    //n_hashes: usize,
    filter: BitVec,
    pub hashers: Vec<NtHasher>, //a vec of hash functions maybe ,or smth like an ntHash build je sais pas
}

impl BloomFilter {
    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u32, k: usize) -> Self {
        Self {
            //size,
            //n_hashes,
            filter: BitVec::from_elem(size, false),
            hashers: init_hashers(n_hashes, seed, k),
        }
    }

    pub fn new(size: usize, n_hashes: usize, k: usize) -> Self {
        let seed: u32 = 42;
        Self::new_with_seed(size, n_hashes, seed, k)
    }

    ///checks if the kmer with specified minimizer hash, and multiple hashes is
    ///inside the bloom filter, inserts it if needed
    pub fn check_and_insert(&mut self, hashed_minimizer: u64, kmer_s_hashes: Vec<u32>) -> bool {
        let mut present: bool = true;
        
        for hash in kmer_s_hashes {
            //to get the address, heavy bits are from the minimizer (giving the block)
            //and light bits are given by the hash of the kmer himself
            let address: usize = 
                (hashed_minimizer as usize)*BLOCK_SIZE + (hash as usize%BLOCK_SIZE);
            if !self.filter.get(address).unwrap_or(false) {
                self.filter.set(address, true);
                present = false;
            }
        }
        present
    }

    pub fn check_true_bits(&self) -> usize {
        let mut counter: usize = 0;
        for i in 0..4294967296 {
            if self.filter.get(i).unwrap() {
                counter += 1;
            }
        }
        counter
    }
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
