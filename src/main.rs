//! kmer counter that uses a blocked bloom filter with hashing of minimizers to determine blocks
//! is tuned for a specific labtop (for now) and only supports up to 31-mers (and not optimal for
//! k<31)

mod input;
mod bloom;
mod counter;
mod output;
mod unit_tests_one_day;
pub mod utils;
pub mod decyclers;
pub mod super_bitvec;
pub mod minimizers;

use input::{read_fof, read_fasta, Hell};
use minimizers::{decycling_mins_x_pos, minimizers_x_positions};
use decyclers::{Decycler};
use bloom::BloomFilter;
use counter::{CountTable};
use utils::{xorshift_u64, xorshift_u128};
use output::{write_output};
use packed_seq::{Seq, PackedSeqVec, SeqVec, PackedSeq};
use std::env; //for backtrace
use rayon::prelude::*;
use clap::Parser;
use std::sync::Mutex;
use needletail::parse_fastx_file;
use rand::Rng;
use std::time::{Duration, Instant};


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

    ///length of the minimizers for grouping the kmers, has to be inferior to k and inferior to 32
    #[arg(short, long, default_value_t = 11)]
    m: u16,

    ///length of the lmers, needs to be inferior or equal to k
    #[arg(short, long, default_value_t = 31)]
    l: u16,

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
    ///warning : deprecated, probably not functionnal anymore, especially for calc false pos rates
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

    ///enables using simd_minimizer instead of any weird thing i cooked up myself
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    simd_minimizer: bool,

}

pub fn main() {
    //for debug
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }
    //defining all variables constants that are based on the argument input
    let args = Args::parse();
    let k: u16 = args.k;
    let m: u16 = args.m;
    let l: u16 = args.l;
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

    assert!(k >= l); //cant have the lmers be longer than the kmers

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

    let filename = args.input;

    //we create the needed data structures to store everything
    let bloom = BloomFilter::new(size, n_hashes, k as usize, block_size, nb_blocks);
    let hash_table = CountTable::new(table_size, table_block_size);

    //calculating decycling sets as well as its time
    let mut duration_overhead_decycling: Duration;
    let mut decycler_set: Decycler;
    let debut = Instant::now();
    if !args.simd_minimizer {
        decycler_set = Decycler::new(m);
        decycler_set.compute_blocks();
    } else {
        decycler_set = Decycler::new(1);
    }
    duration_overhead_decycling = debut.elapsed();

    //anti optims variable
    let kmer_sum: Mutex<u64> = Mutex::new(0);

    //used to check for false negative rate at the end
    let false_neg_list: Mutex<Vec<PackedSeqVec>> = Mutex::new(Vec::new());

    //now we parse and treat each input method
    if args.input_type == 0 {
        let iter_files = read_fof(filename.to_string());
        iter_files.chunks(sequential_fallback).par_bridge().for_each(|chunk| {
            let mut all_addresses: Vec<usize> = vec![0; 7*(2*k-m) as usize];
            for line in chunk {
                let sequence = read_fasta(line.to_string());
                let local_kmer_sum =
                    handle_sequence(&bloom, &hash_table, sequence, k, m, nb_blocks, 
                    one_to_one, no_bloom, no_hashtable, &mut all_addresses, &decycler_set, l);
                let mut total_sum = kmer_sum.lock().unwrap();
                *total_sum = total_sum.wrapping_add(local_kmer_sum);
                drop(total_sum);
            }
        });
    } else if args.input_type == 1 {
        let reader = parse_fastx_file(&filename).expect("valid path/file");

        let chunked_lines = Hell {
            fxreader : reader,
            chunk_size : sequential_fallback,
        };

        chunked_lines.par_bridge().for_each(|chunk| {
            let mut block_lines_counter: usize = 0;

            //to reduce number of created and destroyed vectors throughout
            let mut all_addresses: Vec<usize> = vec![0; 7*(2*k-m) as usize];

            for line in chunk {
                // /!\/!\ assuming single line writing, so that each line corresponds to a
                // sequence
                //with this assumption make a packedseq from the sequence
                let sequence = PackedSeqVec::from_ascii(&line);

                //roll a dice to add to the false negatives checker
                let dice_roll = rand::rng().random_range(0..5000);
                if dice_roll == 0 {
                    let mut false_negs = false_neg_list.lock().unwrap();
                    false_negs.push(sequence.clone());
                    drop(false_negs);
                }

                if only_parse {
                    block_lines_counter += sequence.len();
                } else if sequence.len() >= k as usize+2 {
                    let local_kmer_sum =
                        handle_sequence(&bloom, &hash_table, sequence, k, m, nb_blocks,
                        one_to_one, no_bloom, no_hashtable, &mut all_addresses, &decycler_set, l);
                    if no_bloom {
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
    let false_negs = false_neg_list.lock().unwrap().to_vec();
    let (false_negative_rate, false_positive_rate) = bloom.count_false_bloom(false_negs, k, m, l, &decycler_set);
    if !args.auto_bench {
        println!("false negative rate : {false_negative_rate}");
        println!("false positive rate : {false_positive_rate}");
    }
    
    //to prevent optims
    if only_parse {
        let total_counter = kmer_sum.lock().unwrap();
        print!("Antim optim counter {total_counter}");
    }



    //maintenant on s'occupe de la sortie et tout la
    let final_count: Vec<u64> = hash_table.calculate_output();
    let _ = write_output(&final_count);


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
            false_negative_rate,
            false_positive_rate,
            duration_overhead_decycling,
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
                let (n_z_bloom, max_bloom, median_bloom, average_bloom, fill_counter) = bloom.count_it_all();
                let n_z_bloom_rate: f64 = n_z_bloom as f64/nb_blocks as f64;
                let max_bloom_rate: f64 = max_bloom as f64/block_size as f64;
                let median_bloom_rate: f64 = median_bloom as f64/block_size as f64;
                let average_bloom_rate: f64 = average_bloom as f64/block_size as f64;
                let overfilled_rate: f64 = fill_counter as f64/n_z_bloom as f64;

                println!("Non zero bf amount : {n_z_bloom}");
                println!("Non zero bloom filter block rates : {n_z_bloom_rate}");
                println!("Max bloom fill rate : {max_bloom_rate}");
                println!("Median fill rate : {median_bloom_rate}");
                println!("Average fill rate : {average_bloom_rate}");
                println!("Overfilled rate : {overfilled_rate}");
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
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    nb_blocks: usize,
    one_to_one: bool,
    no_bloom: bool,
    no_hashtable: bool,
    all_addresses: &mut Vec<usize>,
    decycler_set: &Decycler,
    l: u16,
    ) -> u64 {
    if original_sequence.len() <= k as usize+2 {
        return 0;
    }
    let address_mask: usize = bloom.block_size-1;
    let mut local_kmer_sum: u64 = 0;
    let (super_kmers_positions, minimizer_values, sequence): (Vec<u32>, Vec<u64>, PackedSeqVec);
    if decycler_set.m > 1 {
        //we use the decycler
        (super_kmers_positions, minimizer_values, sequence) =
            decycling_mins_x_pos(original_sequence, k, m, decycler_set);
    } else {
        //we use simd_minimizer
        (super_kmers_positions, minimizer_values, sequence) =
            minimizers_x_positions(original_sequence, k, m);
    }

    //quick check that we don't have abherrent results
    assert!(super_kmers_positions.len()==minimizer_values.len(), 
        "Superkmers and minimizers have different length.");

    let mut kmer_number: usize = 0;
    //compute all hashes at once to g faster than computing them 1 by 1
    for i in 0..super_kmers_positions.len()-1 {
        //using minimizer hashing for now to be sure its not a source of problems, will see if
        //removing it doesn't break anything later
        let hashed_minimizer: u64;
        if one_to_one {
            hashed_minimizer = minimizer_values[i]&(nb_blocks as u64-1);
        } else {
            hashed_minimizer = xorshift_u64(minimizer_values[i])&(nb_blocks as u64-1);
        }
        if no_bloom {
            //prevent optims
            local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
        } else {
            if l <= 31 {
                kmer_number = 
                    handle_super_kmer(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                    bloom, hash_table, k, hashed_minimizer, kmer_number,
                    no_hashtable, all_addresses, address_mask, l);
            } else {
                kmer_number = 
                    handle_super_kmer_u128(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                    bloom, hash_table, k, hashed_minimizer, kmer_number,
                    no_hashtable, all_addresses, address_mask, l);
            }
        }
    }
    //not forgetting the last element of the list
    let hashed_minimizer: u64;
    if one_to_one {
        hashed_minimizer = minimizer_values[minimizer_values.len()-1]&(nb_blocks as u64-1);
    } else {
        hashed_minimizer = xorshift_u64(minimizer_values[minimizer_values.len()-1])&(nb_blocks as u64-1);
    }
    if no_bloom {
        //prevent optims
        local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
    } else {
        if l <= 31 {
            let _ = 
                handle_super_kmer(super_kmers_positions[super_kmers_positions.len()-1], 
                (sequence.len()+1-k as usize) as u32,
                &sequence, 
                bloom, hash_table, k, hashed_minimizer,
                kmer_number, no_hashtable, all_addresses, address_mask, l);
        } else {
            let _ = 
                handle_super_kmer_u128(super_kmers_positions[super_kmers_positions.len()-1], 
                (sequence.len()+1-k as usize) as u32,
                &sequence,
                bloom, hash_table, k, hashed_minimizer,
                kmer_number, no_hashtable, all_addresses, address_mask, l);
        }
    }

    //is here only to prevent optimisations in case no bloom filters
    local_kmer_sum
}

fn handle_super_kmer(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec,
    bloom: &BloomFilter, 
    _hash_table: &CountTable, 
    k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, _no_hashtable: bool, all_addresses: &mut Vec<usize>, 
    address_mask: usize, l: u16) -> usize {
    let mut last_relevant_index: usize = 0;
    for j in (start_pos as usize)..(end_pos as usize) + (k-l) as usize{
        let lmer: PackedSeq = sequence.slice(j..j+l as usize);
        let mut hash: u64 = xorshift_u64(lmer.as_u64());


        for _i in 0..bloom.n_hashes {
            let address = hash as usize & address_mask;
            all_addresses[last_relevant_index] = address;
            last_relevant_index += 1;
            hash = xorshift_u64(hash);
        }

        kmer_number+=1;

    }
    let relevant_addresses = &mut all_addresses[..last_relevant_index];
    let blocknum: usize = (hashed_minimizer as usize)&1023;
    let subblocknum: usize = ((hashed_minimizer as usize)>>10)&((bloom.nb_blocks>>10)-1);
    let mut block = bloom.filter[blocknum].lock().unwrap();
    let subblock = &mut block[subblocknum];
    for address in relevant_addresses {
            if !subblock.get(*address) {
                subblock.set(*address, true);
            }
    }
    drop(block);
    kmer_number
}

fn handle_super_kmer_u128(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec,
    bloom: &BloomFilter, 
    _hash_table: &CountTable, 
    k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, _no_hashtable: bool, all_addresses: &mut Vec<usize>, 
    address_mask: usize, l: u16) -> usize {
    let mut last_relevant_index: usize = 0;
    for j in (start_pos as usize)..(end_pos as usize) + (k-l) as usize {
        let lmer: PackedSeq = sequence.slice(j..j+l as usize);
        let mut hash: u128 = xorshift_u128(lmer.as_u128());


        for _i in 0..bloom.n_hashes {
            let address = hash as usize & address_mask;
            all_addresses[last_relevant_index] = address;
            last_relevant_index += 1;
            hash = xorshift_u128(hash);
        }

        kmer_number+=1;

    }
    let relevant_addresses = &mut all_addresses[..last_relevant_index];
    let blocknum: usize = (hashed_minimizer as usize)&1023;
    let subblocknum: usize = ((hashed_minimizer as usize)>>10)&((bloom.nb_blocks>>10)-1);
    let mut block = bloom.filter[blocknum].lock().unwrap();
    let subblock = &mut block[subblocknum];
    for address in relevant_addresses {
            if !subblock.get(*address) {
                subblock.set(*address, true);
            }
    }
    drop(block);
    kmer_number
}


fn write_auto_bench_stdout(
    no_bloom : bool, 
    no_hashtable : bool,
    bloom: BloomFilter,
    hash_table: CountTable,
    nb_blocks: usize,
    block_size: usize,
    table_size: usize,
    table_block_size: usize,
    false_positive_rate: f64,
    false_negative_rate: f64,
    duration_overhead_decycling: Duration,
    ) {
    let mut print_string = String::new();
    //writes every number looked for by the benchmark programm in a single line
    //also does all the counting
    if !no_bloom {
        let (n_z_bloom, max_bloom, median_bloom, average_bloom, fill_counter) = bloom.count_it_all();
        let n_z_bloom_rate: f64 = n_z_bloom as f64/nb_blocks as f64;
        let max_bloom_rate: f64 = max_bloom as f64/block_size as f64;
        let median_bloom_rate: f64 = median_bloom as f64/block_size as f64;
        let average_bloom_rate: f64 = average_bloom as f64/block_size as f64;
        let overfilled_rate: f64 = fill_counter as f64/n_z_bloom as f64;

        print_string += 
            &format!("{n_z_bloom_rate}|{max_bloom_rate}|{average_bloom_rate}|{median_bloom_rate}|{overfilled_rate}");
    } else {
        print_string += &format!("0|0|0|0");
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

    //false negatives and false potitives rates
    print_string += &format!("|{:.3}|{:.3}", false_positive_rate, false_negative_rate);

    //duration of overhead decycling set calculation
    print_string += &format!("|{}", duration_overhead_decycling.as_secs());

    println!("{print_string}");
}
