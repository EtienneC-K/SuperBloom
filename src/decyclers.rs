//! module to implement decycling minimizers
//! because we will always use (relatively) small minimizers
//! the strategy will always be to first compute the complete list of all decycled minimizers, and
//! then to query it in read only accross the different threads

use rayon::prelude::*;
//use packed_seq::{PackedSeqVec, SeqVec, PackedSeq, Seq};
use packed_seq::{PackedSeq, Seq};

const CYCLER_BLOCK_SIZE: usize = 512; //just always keep it a power of 2

pub struct Decycler {
    m: u16, //minimizer length
    pub direct_list: Vec<Vec<u64>>, //stores in "booleans" if a string is a direct decycling minimizer
                                //based on its address
    //indirect_lit, to add
}

impl Decycler {

    ///initialization, creating enough space for all the minimizers
    pub fn new(m: u16) -> Self {
        let direct_list: Vec<Vec<u64>> = vec![vec![0; CYCLER_BLOCK_SIZE]; (1<<(2*m))/(64*CYCLER_BLOCK_SIZE)];
        init_vec_ci(m);
        Self {
            m,
            direct_list,
        }
    }

    ///computes the belonging (or not) of all the kmers
    pub fn compute_blocks(&mut self) {
        let vec_ci: Vec<f64> = init_vec_ci(self.m);
        self.direct_list.par_iter_mut().enumerate().for_each(|(i, mut block)| {
            //println!("ca compute un block en léééégende");
            compute_block(i, &mut block, self.m, &vec_ci);
        })
    }

    pub fn lookup(&self, minimizer: PackedSeq) -> bool {
        //start by converting the kmer to an address we can use
        let address = minimizer.as_u64() as usize;

        //finnding the block
        let block_adress: usize = address/(CYCLER_BLOCK_SIZE*64);
        let block = &self.direct_list[block_adress];

        //lookup the corresponding u64
        let integer_adress: usize = (address%(CYCLER_BLOCK_SIZE*64))/64;
        let integer: u64 = block[integer_adress];

        //reading the correct bit using bitshifting
        let boolean: bool = if (integer>>(63-address%64))%2 == 1 {true} else {false};

        //return
        boolean
    }
}


fn compute_block(i: usize, block: &mut Vec<u64>, m: u16, vec_ci: &Vec<f64>) {
    let mut kmer: u64 = (i*CYCLER_BLOCK_SIZE*64) as u64;

    for j in 0..CYCLER_BLOCK_SIZE {
        let mut to_insert: u64 = 0;
        for k in 0..64 {
            let is_decycler: bool = compute_membership(kmer, m, vec_ci);
            if is_decycler {
                to_insert += 1<<(63-k);
            }

            kmer +=1;
        }
        block[j] = to_insert;
    }
}

///algorithm to check if a kmer is member of a minimum decysling set
///see "Efficient minimizer orders for large values of k using minimum decycling sets" by
///David Pellow, Lianrong Pu, Baris Ekim, Kior Kotlar, Bonnie Berger, Ron Shamir and Yaron
///Orenstein; page 3
pub fn compute_membership(kmer: u64, m: u16, vec_ci: &Vec<f64>) -> bool {
    let epsilon: f64 = 0.00001;
    let mut imaginary_x: f64 = 0.0; 
    let mut imaginary_x_prime: f64 = 0.0;
    for i in 0..m {
        let x_i: u64 = (kmer>>(2*(m-i-1)))%4;
        imaginary_x += vec_ci[i as usize]*x_i as f64;
        let i_prime: usize = if i<m-1 {i as usize+1} else {0};
        imaginary_x_prime += vec_ci[i_prime]*x_i as f64;
    }
    //println!("partie imaginaire : {imaginary_x}, du précédent : {imaginary_x_prime}");

    if imaginary_x > epsilon {
        if imaginary_x_prime <= epsilon {
            return true
        }
    }else if imaginary_x >= -epsilon && imaginary_x <= epsilon { //testing equality actually
        if imaginary_x_prime >= -epsilon && imaginary_x_prime <= epsilon {
            let mut k: u16 = 0;
            for l in 0..2*m {
                let x_l_mod_m: u64 = (kmer>>(2*(m-(l%m)-1)))%4;
                let x_k: u64 = (kmer>>(2*(m-k-1)))%4;
                if x_l_mod_m < x_k {return false};
                if x_l_mod_m > x_k {k = 0} else {k += 1};
                if (l >= m-1) && (k%m == 0) {return true};
            }
        }
    }
    return false
}

pub fn init_vec_ci(m: u16) -> Vec<f64> {
    let mut vec_ci: Vec<f64> = Vec::with_capacity(m as usize);
    for i in 0..m {
        vec_ci.push((2.0*std::f64::consts::PI*i as f64/m as f64).sin());
    }
    vec_ci
}
