///module that contains all necessary functions for the actual counting part of the kmers once they
///have passed through the bloom filter

//use std::ptr;
use bitvec::prelude::*;
//use packed_seq::packed_seq::PackedSeqBase;
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};

///hash table taht will store all the kmers (but not yet their count)
pub struct CountTable {
    table: Vec<u64>, //ca store up to 31-mers bc of the chosen size
    counters: Vec<u32>,
    skip_counter: u64, //counts the amount of kmers that were not inserted
    //zzero: PackedSeqVec,
}

impl CountTable {
    //const TABLE_SIZE: usize = 450000000; //450 millions
    const TABLE_SIZE: usize = 3200000; //3.2 millions
    const MAX_RETRIES: usize = 10;
    
    pub fn new() -> Self {
        let table: Vec<u64> = vec![u64::MAX; Self::TABLE_SIZE];
        let counters: Vec<u32> = vec![0; Self::TABLE_SIZE];
        let skip_counter: u64 = 0;
        Self {
            table,
            counters,
            skip_counter,
        }
    }

    ///checks if the kmer is already inserted, or inserts it if its not, and then increments
    ///its counter, if after max_retries there is still no place that was found for the kmer
    ///we increment the skip_counter instead
    pub fn insert(&mut self, kmer: u64, hashed_kmer: u32) {
        let mut inserted: bool = false;
        let mut i: usize = 0;
        while i<Self::MAX_RETRIES && !inserted {
            let current_address = ((hashed_kmer as usize) + (i+i.pow(2))/2) % Self::TABLE_SIZE;
            if self.table[current_address] == kmer {
                self.counters[current_address] = self.counters[current_address].saturating_add(1);
                inserted = true;
            }
            //we check the last bit that corresponds to insertion or not
            else if self.table[current_address] == u64::MAX { //checking if unused
                self.table[current_address] = kmer;
                self.counters[current_address] = self.counters[current_address].saturating_add(1);
                inserted = true;
            }
            //lets not forget to increment i
            i+=1;
        }
        if !inserted {
            //that means we skip it
            self.skip_counter = self.skip_counter.saturating_add(1);
        }
    }

    ///function that counts the number of each occurence on a 256 long u64 vector, any occurence
    ///higher than 255 being stored a the last cell
    pub fn calculate_output(&self) -> Vec<u64> {
        let mut final_count_vec: Vec<u64> = vec![0; 256];
        for count in &self.counters {
            if *count > 255 {
                final_count_vec[255] = final_count_vec[255].saturating_add(1);
            } else {
                final_count_vec[*count as usize] = 
                    final_count_vec[*count as usize].saturating_add(1);
            }
        }
        final_count_vec
    }


    //fonction qui sert pour réaliser "l'histogramme" de sorti
}
