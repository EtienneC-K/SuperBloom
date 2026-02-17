///module taht implements my super u64 bitvecs
use std::fmt::Error;
use std::fmt::Display;
use std::fmt;

pub struct SuperBitVec {
    vector: Vec<u64>,
    size: usize,
}

impl SuperBitVec {

    ///creates a new superbitvec, filled with 0's, "makes the actual complete structure"
    pub fn new(size : usize) -> Self {
        let vec_size: usize = if size%64==0 {size/64} else {size/64+1};
        let vector: Vec<u64> = vec![0; vec_size];
        
        Self {
            vector,
            size,
        }
    }

    //sets the bit at "index" address to value
    pub fn set(&mut self, address: usize, value: bool) {
        //first check the address is legal
        if address >= self.size {
            panic!("Index out of SuperBitVec range");
        }

        //compute which bit to insert
        //let block_num: usize = address/64;
        //let mut to_insert: u64 = 1<<(63-address%64) as u64; //trust the calculation
        //REMOVED MODULO
        let block_num: usize = address>>6;
        let mut to_insert: u64 = 1<<(63-(address&63)) as u64; //trust the calculation
        if value == false {
            to_insert = u64::MAX-to_insert;
        }

        //now performing the insertion with an atomic or
        if value == true {
            self.vector[block_num] = self.vector[block_num] | to_insert;
        } else {
            self.vector[block_num] = self.vector[block_num] & to_insert;
        }
    }

    ///getter for a certain bit
    pub fn get(&self, address: usize) -> bool {
        //let block = self.vector[address/64];
        //let boolean: bool = if (block>>(63-address%64))%2 == 1 {true} else {false};
        //REMOVED MODULO
        let block = self.vector[address>>6];
        let boolean: bool = if (block>>(63-(address&63)))&1 == 1 {true} else {false};
        boolean
    }

    ///len method cuz its used quite a bit
    pub fn len(&self) -> usize {
        self.size
    }
}


impl Display for SuperBitVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), Error> {

        //is meant for some debugging on smaller examples, so we have a max threshhold
        let write_threshhold: usize = 8192;
        if self.size > write_threshhold {
            return write!(f, "SuperBitVec of length over {write_threshhold} wasn't written.");
        }

        //special case for the first block depends on the actual length of the SuperBitVec
        let mut to_write = String::new();

        //let to_push: usize = if self.size%64==0 {0} else {64-self.size%64};
        //REMOVED MODULO
        let to_push: usize = if self.size&63==0 {0} else {64-(self.size&63)};
        let first_block = self.vector[0] >> to_push;
        to_write += &format!("{:0width$b}", first_block, width = 64-to_push); 

        //let last_block_number: usize = if self.size%64==0 {self.size/64} else {self.size/64+1};
        //REMOVED MODULO
        let last_block_number: usize = if self.size&63==0 {self.size>>6} else {(self.size>>6)+1};
        for i in 1..(last_block_number) {
            to_write += &format!("{:064b}", self.vector[i]);
        }
        write!(f, "{to_write}")
    }
}

impl Clone for SuperBitVec {
    fn clone(&self) -> Self {
        let vector = self.vector.clone();
        let size = self.size;
        Self {
            vector,
            size,
        }
    }
}
