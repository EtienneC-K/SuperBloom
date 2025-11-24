//! kmer counter that uses a blocked bloom filter with hashing of minimizers to determine blocks
//! is tuned for a specific labtop (for now) and only supports up to 31-mers (and not optimal for
//! k<31)

mod input;
mod minimizers;
mod bloom;
mod utils;
mod counter;
mod output;

use input::{read_fof, read_fasta};
use minimizers::minimizers_x_positions;
use bloom::{BloomFilter, BLOCK_SIZE, NB_BLOCKS};
use counter::{CountTable};
use utils::{xorshift_u64};
use output::write_output;
use seq_hash::{KmerHasher};
use packed_seq::{Seq, PackedSeqVec, SeqVec, PackedSeq};
use bitvec::prelude::*;
use std::env; //for backtrace

fn main() {
    //for debug
    unsafe {
        env::set_var("RUST_BACKTRACE", "1");
    }
    //read arguments
    //TODO take arguments instead of having constants over here
    let k: u16 = 31;
    let m: u16 = 11;
    let n_hashes: usize = 7; //static 7 for now to aim for <1% false positives
    //let nb_blocks: usize = 1<<14; //16 384 for now, will see later to make it varaible
    //let size: usize = 1<<35; // 34 359 738 368bits so 4 294 967 296bytes
    //let filename = read_arguments
    let size: usize = 1<<27; // 34 359 738 368bits so 4 294 967 296bytes
    //let nb_blocks: usize = 1<<15; //16 384 for now, will see later to make it varaible

    //for now check of size awith the blocks, later only two of them will be specified
    assert!(size == BLOCK_SIZE*NB_BLOCKS, "Error on filter and block sizes, do not match.");

    let filename = "fasta_reads/listing.txt";

    //we create the needed data structures to store everything
    let mut bloom = BloomFilter::new(size, n_hashes, k as usize);
    let mut hash_table = CountTable::new();

    let iter_files = read_fof(filename.to_string());
    //pin_mut!(iter_files); //needed for iteration

    for fasta_name in iter_files {
        let sequence = read_fasta(fasta_name);
        handle_fasta(&mut bloom, &mut hash_table, sequence, k, m, n_hashes, NB_BLOCKS);
    }


    //maintenant on s'occupe de la sortie et tout la
    let final_count: Vec<u64> = hash_table.calculate_output();
    let _ = write_output(&final_count);

}

fn handle_fasta(
    bloom: &mut BloomFilter,
    hash_table: &mut CountTable,
    sequence: PackedSeqVec,
    k: u16,
    m: u16,
    n_hashes: usize,
    nb_blocks: usize,
    ) {
    let (super_kmers_positions, minimizer_values, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec)
                                                    = minimizers_x_positions(sequence, k, m);

    //super_kmers_positions, minimizer_values, sequence = minimizers_x_positions(sequence, k, m);
    //quick check that we don't have abherrent results
    assert!(super_kmers_positions.len()==minimizer_values.len(), 
        "Superkmers and minimizers have different length.");



    //hasher for the minimizers, will probably be removed later on

    let mut kmer_number: usize = 0;
    //compute all hashes at once to g faster than computing them 1 by 1
    let mut all_hashes: Vec<Vec<u32>> = Vec::new();
    for i in 0..bloom.hashers.len() {
        all_hashes.push(bloom.hashers[i].hash_kmers_simd(sequence.as_slice(), 1).collect());
    }
    for i in 0..super_kmers_positions.len()-1 {
        //using minimizer hashing for now to be sure its not a source of problems, will see if
        //removing it doesn't break anything later
        let hashed_minimizer = xorshift_u64(minimizer_values[i])%(nb_blocks as u64);
        //cf magnifique dessin de quels kmers appartienent à quel super_kmer
        kmer_number = 
            handle_super_kmer(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
            n_hashes, bloom, hash_table, k, hashed_minimizer, &all_hashes, kmer_number);
    }
    //pas oublier le dernier morceau de la liste a évaluer maintenant
    let hashed_minimizer = 
        xorshift_u64(minimizer_values[minimizer_values.len()-1])%(nb_blocks as u64);
    let _ = 
        handle_super_kmer(super_kmers_positions[super_kmers_positions.len()-1], 
        (sequence.len()-1-k as usize) as u32,
        &sequence, 
        n_hashes, bloom, hash_table, k, hashed_minimizer,
        &all_hashes, kmer_number);
}

fn handle_super_kmer(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec, n_hashes: usize,
    bloom: &mut BloomFilter, hash_table: &mut CountTable, k: u16, hashed_minimizer: u64, 
    all_hashes: &Vec<Vec<u32>>, mut kmer_number: usize) -> usize {
    for j in (start_pos as usize)..(end_pos as usize) {
        let kmer: PackedSeq = sequence.slice(j..j+k as usize);
        let mut kmer_s_hashes: Vec<u32> = Vec::new();

        for i2 in 0..n_hashes {
            kmer_s_hashes.push(all_hashes[i2][kmer_number]);
        }

        let already_in = bloom.check_and_insert(hashed_minimizer, kmer_s_hashes);
        //do_smth if it was already in, like adding it to a hash_table for counting
        //problem with that : its gonna take an awful lot of space i think (it does)
        if already_in {
            let kmer_hash = all_hashes[0][kmer_number];
            //let bitvec_kmer: BitVec = convert_seqkmer(kmer);
            hash_table.insert(kmer.as_u64(), kmer_hash); //we take the first hash for the hash
                                                            //table as well
        }
        kmer_number+=1;
    }
    kmer_number
}
