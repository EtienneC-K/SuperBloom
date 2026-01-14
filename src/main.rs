//! kmer counter that uses a blocked bloom filter with hashing of minimizers to determine blocks
//! is tuned for a specific labtop (for now) and only supports up to 31-mers (and not optimal for
//! k<31)

mod input;
mod minimizers;
mod bloom;
mod utils;
mod counter;
mod output;

use input::{read_fof, read_fasta, read_lines, Hell};
use minimizers::minimizers_x_positions;
//use bloom::{BloomFilter, BLOCK_SIZE, NB_BLOCKS};
use bloom::BloomFilter;
use counter::{CountTable};
use utils::{xorshift_u64, xorshift_u32};
use output::{write_output};
use seq_hash::{KmerHasher};
use packed_seq::{Seq, PackedSeqVec, SeqVec, PackedSeq};
use std::env; //for backtrace
use rayon::prelude::*;
use clap::Parser;
use itertools::Itertools;
use std::io::BufReader;
use std::fs::File;
use std::io::BufRead;
use std::sync::Mutex;
use needletail::parse_fastx_file;

///taking care of all the needed command line arguments
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    ///path to the file of file
    input: String,

    ///output file (defautl is out.csv)
    #[arg(short, long, default_value_t = String::from("out.csv"))]
    output: String,

    ///length of the kmers to count
    #[arg(short, long, default_value_t = 31)]
    k: u16,

    ///length of the minimizers for grouping the kmers, has to be inferior to k
    #[arg(short, long, default_value_t = 11)]
    m: u16,

    ///number of hashes for the bloom filter
    #[arg(short, long, default_value_t = 7)]
    n_hashes: usize,

    ///size (in bits) of the bloom filter, expressed as a power of 2
    #[arg(short, long, default_value_t = 33)]
    size: usize,

    ///size (in bits) of each of the bloom filters blocks, expressed as a power of 2
    #[arg(short, long, default_value_t = 14)]
    block_size: usize,

    ///size (in number of slots) of the hash table for the final count, expressed as a power of 2 
    #[arg(long, default_value_t = 28)]
    table_size: usize,

    ///size (in number of slots) of each of the hash table's blocks, as a power of 2
    #[arg(long, default_value_t = 14)]
    table_block_size: usize,

    ///number of threads (default 1)
    #[arg(short, long, default_value_t = String::from("1"))]
    threads: String,

    ///input method, can be 0 for a file of file directing to single fastas, or 1 for a single
    ///multi_fasta (in which case we assume very hard that lines' lenghts do not cap a 80, but at the 
    ///reads' lengths)
    ///dev note: input_type 0 is way behind in commits, do not use
    #[arg(long, default_value_t = 0)]
    input_type: u8,

    ///argument to set the blocks to match one to one the minimizers, based on minimizer length,
    ///and bloom size and overrides the blocksize
    #[arg(short, long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    one_to_one: bool,

    ///number of reads to be distributed in a row to each thread
    #[arg(long, default_value_t = 100)]
    sequential_fallback: usize,

    //to enable counting outputs, that take time outside of the actual algorithm
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    counting: bool,

    ///to disable all code referring to the hash_table, for testing without it
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_hashtable: bool,

    ///to disable all code referring to the bloom filter, enables no_hashtable
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_bloom: bool,

    ///to disable all code after the parsing part, for bench purposes
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    only_parse: bool,

    ///to change the standard output to just sequence of numbers to be read by a benchmark programm
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    auto_bench: bool,

}

pub fn main() {
    //for debug
    unsafe {
        env::set_var("RUST_BACKTRACE", "1");
    }
    //defining all variables constants that are based on the argument input
    let args = Args::parse();
    let k: u16 = args.k;
    let m: u16 = args.m;
    let n_hashes: usize = args.n_hashes;
    let mut size: usize = 1<<args.size;
    let mut block_size:usize = 1<<args.block_size;
    let mut nb_blocks: usize = size/block_size;
    let mut table_size: usize = 1<<args.table_size;
    let mut table_block_size: usize = 1<<args.table_block_size;
    let one_to_one: bool = args.one_to_one;
    let sequential_fallback: usize = args.sequential_fallback;
    let only_parse: bool = args.only_parse;
    let no_bloom: bool = args.no_bloom || args.only_parse;
    let no_hashtable: bool = args.no_hashtable || args.no_bloom || args.only_parse;

    //for the special case where i want to map 
    if one_to_one {
        nb_blocks = 1<<(2*m);
        block_size = size/nb_blocks;
    }

    if no_bloom {
        size = 1024;
        block_size = 1;
        nb_blocks = 1024;
    }

    if no_hashtable {
        table_size = 1024;
        table_block_size = 1;
    }

    //number of threads allowed
    unsafe {
        env::set_var("RAYON_NUM_THREADS", args.threads);
    }

    //for now check of size awith the blocks, later only two of them will be specified
    //assert!(size == BLOCK_SIZE*NB_BLOCKS, "Error on filter and block sizes, do not match.");

    let filename = args.input;

    //we create the needed data structures to store everything
    let bloom = BloomFilter::new(size, n_hashes, k as usize, block_size, nb_blocks);
    let hash_table = CountTable::new(table_size, table_block_size);

    //anti optims variable
    let kmer_sum: Mutex<u64> = Mutex::new(0);

    //now we parse and treat each input method
    if args.input_type == 0 {
        let iter_files = read_fof(filename.to_string());
        iter_files.chunks(sequential_fallback).par_bridge().for_each(|chunk| {
            for line in chunk {
                let sequence = read_fasta(line.to_string());
                let local_kmer_sum =
                    handle_sequence(&bloom, &hash_table, sequence, k, m, n_hashes, nb_blocks, 
                    one_to_one, no_bloom, no_hashtable);
                let mut total_sum = kmer_sum.lock().unwrap();
                *total_sum = total_sum.wrapping_add(local_kmer_sum);
                drop(total_sum);
            }
        });
    } else if args.input_type == 1 {
        //let file = File::open(filename).unwrap();
        //let reader = BufReader::new(file);
        //let lines = reader.lines();
        let mut reader = parse_fastx_file(&filename).expect("valid path/file");
    

        let chunked_lines = Hell {
            fxreader : reader,
            chunk_size : sequential_fallback,
        };

        chunked_lines.par_bridge().for_each(|chunk| {
        //reader.chunks(sequential_fallback).par_bridge().for_each(|chunk| {
            let mut block_lines_counter: usize = 0;
            for line in chunk {
                //let line = wrapped_line.expect("problem unwrapping a line from the multi fasta");
                //if line.len() > k as usize 
                //    && line.as_bytes()[0] != b'>' 
                //    && line.as_bytes()[0] != b'@' {
                // /!\/!\ assuming single line writing, so that each line corresponds to a
                // sequence
                //with this assumption make a packedseq from the sequence
                let sequence = PackedSeqVec::from_ascii(&line);

                
                if only_parse {
                    block_lines_counter += sequence.len();
                } else {
                    //println!("{:?}", line);
                    let local_kmer_sum =
                        handle_sequence(&bloom, &hash_table, sequence, k, m, n_hashes, nb_blocks,
                        one_to_one, no_bloom, no_hashtable);
                    if !no_bloom {
                        let mut total_sum = kmer_sum.lock().unwrap();
                        *total_sum = total_sum.wrapping_add(local_kmer_sum);
                        drop(total_sum);
                    }
                }

            }
            if only_parse {
                let mut total_sum = kmer_sum.lock().unwrap();
                *total_sum = total_sum.wrapping_add(block_lines_counter as u64);
                drop(total_sum);
            }
        })

    } else {
        panic!("Unrecognized input type, must be 0 or 1");
    }
    
    //to prevent optims
    if only_parse {
        let total_counter = kmer_sum.lock().unwrap();
        print!("Antim optim counter {total_counter}");
    }


    //for fasta_name in iter_files {
    //    let sequence = read_fasta(fasta_name);
    //    handle_fasta(&mut bloom, &mut hash_table, sequence, k, m, n_hashes, NB_BLOCKS);
    //}

    //maintenant on s'occupe de la sortie et tout la
    let final_count: Vec<u64> = hash_table.calculate_output();
    let _ = write_output(&final_count);

    //check des taux de remplissage
    //let taux_bloom = bloom.check_true_bits();
    //let proportion_bloom: f64 = (taux_bloom as f64)/(size as f64);
    //let taux_ht = hash_table.check_filling();
    //let proportion_ht: f64 = (taux_ht as f64)/(table_size as f64);
    //println!("Remplissage du bloom {taux_bloom} ce qui représente une proportion de {proportion_bloom}");
    //println!("Remplissage de la HT {taux_ht} ce qui représente une proportion de {proportion_ht}");

    //printing only a line for the benchmark evaluating programm if option --auto-bench if on
    if args.auto_bench {
        write_auto_bench_stdout(
            no_bloom, 
            no_hashtable,
            bloom,
            hash_table,
            nb_blocks,
            block_size,
            table_size,
            table_block_size,
            )
    }
    else {
        println!("------------------------------------------------------------");
        println!("");
        println!("Parameters : ");
        println!("k : {k}, m: {m}");
        println!("bf size : {size}, block size {block_size}, nb_blocks {nb_blocks}");
        println!("ht size : {table_size}, block size {table_block_size}, nb_blocks {0}", table_size/table_block_size);
        println!("sequetial fallback : {sequential_fallback}");
 
        println!("");
 
        if args.counting {

            if !no_bloom {
                let (n_z_bloom, max_bloom, median_bloom, average_bloom) = bloom.count_it_all();
                let n_z_bloom_rate: f64 = n_z_bloom as f64/nb_blocks as f64;
                let max_bloom_rate: f64 = max_bloom as f64/block_size as f64;
                let median_bloom_rate: f64 = median_bloom as f64/block_size as f64;
                let average_bloom_rate: f64 = average_bloom as f64/block_size as f64;

                println!("Non zero bf amount : {n_z_bloom}");
                println!("Non zero bloom filter block rates : {n_z_bloom_rate}");
                println!("Max bloom fill rate : {max_bloom_rate}");
                println!("Median fill rate : {median_bloom_rate}");
                println!("Average fill rate : {average_bloom_rate}");
            }

            println!("");

            if !no_hashtable {
                let (n_z_ht, max_ht, median_ht, average_ht) = hash_table.count_it_all();
                let n_z_ht_rate: f64 = n_z_ht as f64/(table_size /table_block_size) as f64;
                let max_ht_rate: f64 = max_ht as f64/table_block_size as f64;
                let median_ht_rate: f64 = median_ht as f64/table_block_size as f64;
                let average_ht_rate: f64 = average_ht as f64/table_block_size as f64;

                println!("Non zero ht amount : {n_z_ht}");
                println!("Non zero ht block rates : {n_z_ht_rate}");
                println!("Max ht fill rate : {max_ht_rate}");
                println!("Median ht fill rate : {median_ht_rate}");
                println!("Average ht fill rate : {average_ht_rate}");
            }

            println!("");
        }

        println!("And with all that we get a skip amount of {0}", 
            *hash_table.skip_counter.lock().unwrap());
        println!("");
        let anti_optim_count = *kmer_sum.lock().unwrap();
        println!("anti optim count {anti_optim_count}");
        println!("------------------------------------------------------------");
    }
}

fn handle_sequence(
    bloom: &BloomFilter,
    hash_table: &CountTable,
    sequence: PackedSeqVec,
    k: u16,
    m: u16,
    n_hashes: usize,
    nb_blocks: usize,
    one_to_one: bool,
    no_bloom: bool,
    no_hashtable: bool,
    ) -> u64 {
    if sequence.len() <= k as usize {
        return 0;
    }
    let mut local_kmer_sum: u64 = 0;
    let (super_kmers_positions, minimizer_values, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec)
                                                    = minimizers_x_positions(sequence, k, m);

    //super_kmers_positions, minimizer_values, sequence = minimizers_x_positions(sequence, k, m);
    //quick check that we don't have abherrent results
    assert!(super_kmers_positions.len()==minimizer_values.len(), 
        "Superkmers and minimizers have different length.");



    //hasher for the minimizers, will probably be removed later on

    let mut kmer_number: usize = 0;
    //compute all hashes at once to g faster than computing them 1 by 1
    for i in 0..super_kmers_positions.len()-1 {
        //using minimizer hashing for now to be sure its not a source of problems, will see if
        //removing it doesn't break anything later
        let hashed_minimizer: u64;
        if one_to_one {
            hashed_minimizer = minimizer_values[i]%(nb_blocks as u64);
        } else {
            hashed_minimizer = xorshift_u64(minimizer_values[i])%(nb_blocks as u64);
        }
        //cf magnifique dessin de quels kmers appartienent à quel super_kmer
        if no_bloom {
            //prevent optims
            local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
        } else {
            kmer_number = 
                handle_super_kmer(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                n_hashes, bloom, hash_table, k, hashed_minimizer, kmer_number,
                no_hashtable);
        }
    }
    //pas oublier le dernier morceau de la liste a évaluer maintenant
    let hashed_minimizer: u64;
    if one_to_one {
        hashed_minimizer = minimizer_values[minimizer_values.len()-1]%(nb_blocks as u64);
    } else {
        hashed_minimizer = xorshift_u64(minimizer_values[minimizer_values.len()-1])%(nb_blocks as u64);
    }
    if no_bloom {
        //prevent optims
        local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
    } else {
        let _ = 
            handle_super_kmer(super_kmers_positions[super_kmers_positions.len()-1], 
            (sequence.len()-1-k as usize) as u32,
            &sequence, 
            n_hashes, bloom, hash_table, k, hashed_minimizer,
            kmer_number, no_hashtable);
    }

    //is here only to prevent optimisations in case no bloom filters
    local_kmer_sum
}

fn handle_super_kmer(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec, n_hashes: usize,
    bloom: &BloomFilter, hash_table: &CountTable, k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, no_hashtable: bool) -> usize {
    for j in (start_pos as usize)..(end_pos as usize) {
        let kmer: PackedSeq = sequence.slice(j..j+k as usize);
        let mut kmer_s_hashes: Vec<u64> = Vec::new();
        let mut last_hash: u64 = xorshift_u64(kmer.as_u64());

        for i in 0..bloom.n_hashes-1 {
            last_hash = xorshift_u64(last_hash);
            kmer_s_hashes.push(last_hash);
        }

        let kmer_hash: u64 = kmer_s_hashes[0]; //in case we need a copy of it for insertion
        let already_in = bloom.check_and_insert(hashed_minimizer, kmer_s_hashes);
        //do_smth if it was already in, like adding it to a hash_table for counting
        //problem with that : its gonna take an awful lot of space i think (it does)
        if already_in && !no_hashtable {
            //let bitvec_kmer: BitVec = convert_seqkmer(kmer);
            //we take the first hash for the hash table as it seem to not rlly matter
            hash_table.insert(kmer.as_u64(), kmer_hash, hashed_minimizer);
        }
        kmer_number+=1;
    }
    kmer_number
}

///output function that requires local modules and is therefore written here
fn write_auto_bench_stdout(
    no_bloom : bool, 
    no_hashtable : bool,
    bloom: BloomFilter,
    hash_table: CountTable,
    nb_blocks: usize,
    block_size: usize,
    table_size: usize,
    table_block_size: usize,
    ) {
    let mut print_string = String::new();
    //writes every number looked for by the benchmark programm in a single line
    //also does all the counting
    if !no_bloom {
        let (n_z_bloom, max_bloom, median_bloom, average_bloom) = bloom.count_it_all();
        let n_z_bloom_rate: f64 = n_z_bloom as f64/nb_blocks as f64;
        let max_bloom_rate: f64 = max_bloom as f64/block_size as f64;
        let median_bloom_rate: f64 = median_bloom as f64/block_size as f64;
        let average_bloom_rate: f64 = average_bloom as f64/block_size as f64;

        print_string = 
            format!("{n_z_bloom_rate}|{max_bloom_rate}|{average_bloom_rate}|{median_bloom_rate}");
    } else {
        print_string = format!("0|0|0|0");
    }

    if !no_hashtable {
        let (n_z_ht, max_ht, median_ht, average_ht) = hash_table.count_it_all();
        let n_z_ht_rate: f64 = n_z_ht as f64/(table_size /table_block_size) as f64;
        let max_ht_rate: f64 = max_ht as f64/table_block_size as f64;
        let median_ht_rate: f64 = median_ht as f64/table_block_size as f64;
        let average_ht_rate: f64 = average_ht as f64/table_block_size as f64;

        print_string +=
            &format!("|{n_z_ht_rate}|{max_ht_rate}|{average_ht_rate}|{median_ht_rate}");
    } else {
        print_string += "|0|0|0|0";
    }

    println!("{print_string}");
}
