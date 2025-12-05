///module that contains all necessary functions for the actual counting part of the kmers once they
///have passed through the bloom filter
///module will be thrown out for ahash to see of its faster indeed

use std::sync::Mutex;
use std::collections::HashMap;
use ahash::{AHasher, RandomState};

///hash table taht will store all the kmers (but not yet their count)
pub struct CountTable {
    tables: Vec<Mutex<HashMap<u64, u8, RandomState>>>,
    //table_size: usize,
    ht_nb_blocks: usize,
}

impl CountTable {
    //const TABLE_SIZE: usize = 1<<28;
    //const MAX_RETRIES: usize = 10;
    //const HT_BLOCK_SIZE: usize = 16384;
    //const HT_NB_BLOCKS: usize = Self::TABLE_SIZE/Self::HT_BLOCK_SIZE;
    
    pub fn new(table_size: usize, table_block_size: usize) -> Self {
        let nb_blocks = table_size/table_block_size;
        let mut tables: Vec<Mutex<HashMap<u64, u8, RandomState>>> = Vec::new();
        for _ in 0..nb_blocks {
            tables.push(Mutex::new(HashMap::default()));
        }
        Self {
            tables,
            //table_size,
            ht_nb_blocks: nb_blocks,
        }
    }

    ///checks if the kmer is already inserted, or inserts it if its not, and then increments
    ///its counter, if after max_retries there is still no place that was found for the kmer
    ///we increment the skip_counter instead
    pub fn insert(&self, kmer: u64, hashed_minimizer: u64) {
        let block_address = hashed_minimizer as usize % self.ht_nb_blocks;
        let mut block = self.tables[block_address].lock().unwrap();
        //ahash_insert(block, kmer);
        match block.get_mut(&kmer) {
            Some(count) => *count = count.saturating_add(1),
            None => {let _ = block.insert(kmer, 1);},
        }
    }

    ///function that counts the number of each occurence on a 256 long u64 vector, any occurence
    ///higher than 255 being stored a the last cell
    pub fn calculate_output(&self) -> Vec<u64> {
        let mut final_count_vec: Vec<u64> = vec![0; 256];
        for i in 0..self.tables.len() {
            let block = self.tables[i].lock().unwrap();
            for value in block.values() {
                if *value >= 255 {
                    final_count_vec[255] = final_count_vec[255].saturating_add(1);
                } else {
                    final_count_vec[*value as usize] = final_count_vec[*value as usize].saturating_add(1);
                }
            }
        }
        final_count_vec
    }


    //fonction qui sert pour réaliser "l'histogramme" de sorti
}

fn ahash_insert(ahash_table: HashMap<u64, u8>, kmer: u64) {
    //match ahash_table.get_mut(*kmer) {
    //    Some(count) => *count = *count.saturating_add(1),
    //    None => ahash_table.insert(kmer, 1),
    //}
}
