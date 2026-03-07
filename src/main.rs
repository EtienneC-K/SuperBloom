//! kmer counter that uses a blocked bloom filter with hashing of minimizers to determine blocks
//! is tuned for a specific labtop (for now) and only supports up to 31-mers (and not optimal for
//! k<31)

mod input;
mod bloom;
mod unit_tests_one_day;
pub mod utils;
pub mod decyclers;
pub mod super_bitvec;
pub mod minimizers;

use input::{Hell};
use minimizers::{decycling_mins_x_pos, minimizers_x_positions};
use decyclers::{Decycler};
use bloom::BloomFilter;
use utils::{xorshift_u64, xorshift_u128, sum_vec_bool};
use packed_seq::{Seq, PackedSeqVec, SeqVec, PackedSeq};
use std::env; //for backtrace
use rayon::prelude::*;
use clap::Parser;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use needletail::parse_fastx_file;
use rand::Rng;
use std::time::{Duration, Instant};


///taking care of all the needed command line arguments, first the more open ones, and then the
///"expert" ones
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    ///path to the file to fill the bloom filter
    input: Option<String>,

    ///path to the indexed file to fill the bloom filter
    #[arg(long)]
    indexed_file: Option<String>,

    ///path to a file containing the sequences to be queried
    #[arg(long, default_value_t = String::from(""))]
    query_file: String,

    ///length of the kmers
    #[arg(short, long, default_value_t = 31)]
    k: u16,

    ///quality versus performance parameter, positive integer, usually between 1 and 3, a higher
    ///value will lead to less false positive but slower execution than higher values
    #[arg(short, long, default_value_t = 2)]
    b: u16,

    ///max amount of RAM (integer in GB) to be used, must be at least 1
    #[arg(long, default_value_t = 1)]
    ram: usize,

    ///number of threads to use (default 1)
    #[arg(short, long, default_value_t = String::from("1"))]
    threads: String,

    ///number of hashes for the bloom filter
    #[arg(short, long, default_value_t = 3)]
    n_hashes: usize,

    ///length of the s-mers, needs to be inferior or equal to k, if left at the default value 0
    ///will be set to k-3
    #[arg(short, long, default_value_t = 0)]
    s: u16,

    ///enables the use of "expert" parameters, if this flag is down all the aftermentionned
    ///parameters will
    ///be chosen automatically, if it is up, you can specifiy the ones you want to fine tune
    ///the bloom filters to your exact usage, does require you to know what you're doing to still
    ///have good performances and results.
    /// \n Will override the b and ram basic parameters
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    exp: bool,

    ///length of the minimizers for grouping the kmers, has to be inferior to k and inferior to 32,
    #[arg(short, long, default_value_t = 11)]
    m: u16,

    ///size (in bits) of the bloom filter, expressed as a power of 2, overrides the max ram
    ///parameter
    #[arg(short, long, default_value_t = 33)]
    size: usize,

    ///size (in bits) of each of the bloom filters blocks, expressed as a power of 2
    #[arg(short, long, default_value_t = 13)]
    block_size: usize,

    ///number of reads to be distributed in a row to each thread
    #[arg(long, default_value_t = 100)]
    sequential_fallback: usize,

    //to enable counting outputs, that take time outside of the actual algorithm, but allow to have
    //more insights on the bloom and the choice of settings
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    counting: bool,

    ///to disable all code referring to the bloom filter
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_bloom: bool,

    ///to disable all code after the parsing part, for bench purposes
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    only_parse: bool,

    ///to change the standard output to just sequence of numbers to be read by a benchmark programm
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    auto_bench: bool,

    ///enables using simd_minimizer double decycling minimizers instead of random
    ///is slower but slightly better in terms of false positives
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_simd_minimizer: bool,

}

pub fn main() {
    //for debug
    unsafe {
        env::set_var("RUST_BACKTRACE", "full");
    }
    let whole_run_start = Instant::now();
    //checking the arguments do make some sense
    let args = Args::parse();
    assert!(args.ram >= 1);
    assert!(args.k >= args.s); //cant have the s-mers be longer than the kmers
    assert!(args.m < 32);
    assert!(args.m > 0);

    //defining all variables constants that are based on the argument input
    let k: u16 = args.k;
    let n_hashes: usize = args.n_hashes;
    let s: u16 = if args.s > 0 {args.s} else {k-3};
    assert!(s > 0);
    assert!(s < 62);

    let m: u16;
    let mut size: usize;
    let mut block_size: usize;
    let mut nb_blocks: usize;
    let sequential_fallback: usize;
    let only_parse: bool;
    let no_bloom: bool;
    let no_simd_minimizer: bool;
    let counting: bool;
    let auto_bench: bool;

    if !args.exp {
        //thats for casual users
        size = 1<<((args.ram as f64).log2().floor() as usize)+30;
        block_size = 1<<13;
        nb_blocks = size/block_size;
        m = (nb_blocks as f64).log2() as u16/2 + args.b;
        sequential_fallback = 100;
        counting = false;
        no_bloom = false;
        only_parse = false;
        auto_bench = false;
        no_simd_minimizer = false;


        //if some expert parameters specified but flag is down, warning
        if args.m != 11
        || args.size != 33
        || args.block_size != 13
        || args.sequential_fallback != sequential_fallback
        || args.only_parse != only_parse
        || args.no_bloom != no_bloom
        || args.no_simd_minimizer != no_simd_minimizer
        || args.counting != counting
        || args.auto_bench != auto_bench {
           print!("/!\\ WARNING : experts parameters where specified but the expert parameter ");
           println!("flag is down, are you sure of what you are doing ? Check --help if necessary."); 
        }
    } else {
        //thats only for the experts
        m = args.m;
        size = 1<<args.size;
        block_size = 1<<args.block_size;
        nb_blocks = size/block_size;
        sequential_fallback = args.sequential_fallback;
        only_parse = args.only_parse;
        no_bloom= args.no_bloom || args.only_parse;
        no_simd_minimizer = args.no_simd_minimizer;
        counting = args.counting;
        auto_bench = args.auto_bench;
    }


    if no_bloom {
        size = 1024;
        block_size = 1;
        nb_blocks = 1024;
    }

    //number of threads allowed
    unsafe {
        env::set_var("RAYON_NUM_THREADS", args.threads);
    }

    let filename = resolve_input_path(args.input.as_deref(), args.indexed_file.as_deref());

    //we create the needed data structures to store everything
    let bloom = BloomFilter::new(size, n_hashes, k as usize, block_size, nb_blocks);

    //calculating decycling sets as well as its time
    let duration_overhead_decycling: Duration;
    let mut decycler_set: Decycler;
    let debut = Instant::now();
    if no_simd_minimizer {
        decycler_set = Decycler::new(m);
        decycler_set.compute_blocks();
    } else {
        decycler_set = Decycler::new(1);
    }
    duration_overhead_decycling = debut.elapsed();

    //anti optims variable
    let kmer_sum = AtomicU64::new(0);
    let inserted_kmer_count = AtomicU64::new(0);
    let queried_kmer_count = AtomicU64::new(0);

    //used to check for false negative rate at the end
    let false_neg_list: Mutex<Vec<PackedSeqVec>> = Mutex::new(Vec::new());

    let indexing_start = Instant::now();
    {
        let reader = parse_fastx_file(&filename).expect("valid path/file");

        let chunked_lines = Hell {
            fxreader : reader,
            chunk_size : sequential_fallback,
        };

        chunked_lines.par_bridge().for_each(|chunk| {
            let mut block_lines_counter: usize = 0;
            let mut local_inserted_kmer_count: u64 = 0;
            let mut local_kmer_sum_total: u64 = 0;

            //to reduce number of created and destroyed vectors throughout
            let mut all_addresses: Vec<usize> = vec![0; 7*(2*k-m) as usize];

            for line in chunk {
                // /!\/!\ assuming single line writing, so that each line corresponds to a
                // sequence
                //with this assumption make a packedseq from the sequence
                let sequence = PackedSeqVec::from_ascii(&line);

                //roll a dice to add to the false negatives checker
                if counting {
                    let dice_roll = rand::rng().random_range(0..5000);
                    if dice_roll == 0 {
                        let mut false_negs = false_neg_list.lock().unwrap();
                        false_negs.push(sequence.clone());
                        drop(false_negs);
                    }
                }

                if only_parse {
                    block_lines_counter += sequence.len();
                } else if sequence.len() >= k as usize {
                    local_inserted_kmer_count += (sequence.len() + 1 - k as usize) as u64;
                    let local_kmer_sum =
                        handle_sequence(&bloom, sequence, k, m, nb_blocks,
                        no_bloom, &mut all_addresses, &decycler_set, s);
                    if no_bloom {
                        local_kmer_sum_total = local_kmer_sum_total.wrapping_add(local_kmer_sum);
                    }
                }

            }
            inserted_kmer_count.fetch_add(local_inserted_kmer_count, Ordering::Relaxed);
            if no_bloom {
                kmer_sum.fetch_add(local_kmer_sum_total, Ordering::Relaxed);
            }
            if only_parse {
                kmer_sum.fetch_add(block_lines_counter as u64, Ordering::Relaxed);
            }
        })
    }
    let indexing_duration = indexing_start.elapsed();

    let query_start = Instant::now();
    if args.query_file != "" {
        //this means that we do have to query
        let query_counter: Mutex<usize> = Mutex::new(0);
        let positive_query_counter: Mutex<usize> = Mutex::new(0);

        let reader = parse_fastx_file(&args.query_file).expect("valid path/file");

        let chunked_lines = Hell {
            fxreader : reader,
            chunk_size : sequential_fallback,
        };

        chunked_lines.par_bridge().for_each(|chunk| {
            let mut local_count: usize = 0;
            let mut local_pos_count: usize = 0;
            for line in chunk {
                let sequence = PackedSeqVec::from_ascii(&line);
                let presence_vec: Vec<bool>;
                if s<=31 {
                    presence_vec = bloom.check_sequence(sequence, k, m, s, &decycler_set);
                } else {
                    presence_vec = bloom.check_sequence_u128(sequence, k, m, s, &decycler_set);
                }
                local_count += presence_vec.len();
                local_pos_count += sum_vec_bool(&presence_vec);
            }
            let mut q_count = query_counter.lock().unwrap();
            let mut q_count_pos = positive_query_counter.lock().unwrap();
            *q_count += local_count;
            *q_count_pos += local_pos_count;
            drop(q_count);
            drop(q_count_pos);
            queried_kmer_count.fetch_add(local_count as u64, Ordering::Relaxed);
        });

        if !auto_bench {
            let q_count = query_counter.lock().unwrap();
            let q_count_pos = positive_query_counter.lock().unwrap();
            println!("Number of kmer queried : {q_count}");
            println!("Number of positives : {q_count_pos}");
        }
    }
    let query_duration = query_start.elapsed();

    let (false_negative_rate, false_positive_rate) = if counting {
        let false_negs = false_neg_list.lock().unwrap().to_vec();
        bloom.count_false_bloom(false_negs, k, m, s, &decycler_set)
    } else {
        (0.0, 0.0)
    };
    let whole_run_duration = whole_run_start.elapsed();
    if !auto_bench {
        let inserted_kmers = inserted_kmer_count.load(Ordering::Relaxed);
        let queried_kmers = queried_kmer_count.load(Ordering::Relaxed);
        let whole_run_ns = whole_run_duration.as_nanos() as f64;
        println!("inserted kmers : {inserted_kmers}");
        if inserted_kmers > 0 {
            println!("ns per inserted kmer : {}", whole_run_ns / inserted_kmers as f64);
        } else {
            println!("ns per inserted kmer : N/A");
        }
        if queried_kmers > 0 {
            println!("ns per queried kmer : {}", whole_run_ns / queried_kmers as f64);
        } else {
            println!("ns per queried kmer : N/A");
        }
        println!("total indexing time (s) : {}", indexing_duration.as_secs_f64());
        println!("total query time (s) : {}", query_duration.as_secs_f64());
    }
    
    //to prevent optims
    //printing only a line for the benchmark evaluating programm if option --auto-bench if on
    if auto_bench {
        write_auto_bench_stdout(
            no_bloom, 
            bloom,
            nb_blocks,
            block_size,
            false_negative_rate,
            false_positive_rate,
            duration_overhead_decycling,
            )
    }
    else {
        println!("Parameters : ");
        println!("k : {k}, m: {m}");
        println!("bf size : {size}, block size {block_size}, nb_blocks {nb_blocks}");

        if counting {

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
        }
    }
}

fn resolve_input_path(input: Option<&str>, indexed_file: Option<&str>) -> String {
    match (input, indexed_file) {
        (Some(path), None) => path.to_owned(),
        (None, Some(path)) => path.to_owned(),
        (None, None) => panic!("an input file path is required"),
        (Some(_), Some(_)) => panic!("use either the positional input or --indexed-file, not both"),
    }
}

fn handle_sequence(
    bloom: &BloomFilter,
    original_sequence: PackedSeqVec,
    k: u16,
    m: u16,
    nb_blocks: usize,
    no_bloom: bool,
    all_addresses: &mut Vec<usize>,
    decycler_set: &Decycler,
    l: u16,
    ) -> u64 {
    if original_sequence.len() < k as usize {
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
        let hashed_minimizer: u64 = xorshift_u64(minimizer_values[i])&(nb_blocks as u64-1);
        if no_bloom {
            //prevent optims
            local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
        } else {
            if l <= 31 {
                kmer_number = 
                    handle_super_kmer(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                    bloom, k, hashed_minimizer, kmer_number,
                    all_addresses, address_mask, l);
            } else {
                kmer_number = 
                    handle_super_kmer_u128(super_kmers_positions[i], super_kmers_positions[i+1], &sequence, 
                    bloom, k, hashed_minimizer, kmer_number,
                    all_addresses, address_mask, l);
            }
        }
    }
    //not forgetting the last element of the list
    let hashed_minimizer: u64 = xorshift_u64(minimizer_values[minimizer_values.len()-1])&(nb_blocks as u64-1);
    if no_bloom {
        //prevent optims
        local_kmer_sum = local_kmer_sum.wrapping_add(hashed_minimizer);
    } else {
        if l <= 31 {
            let _ = 
                handle_super_kmer(super_kmers_positions[super_kmers_positions.len()-1], 
                (sequence.len()+1-k as usize) as u32,
                &sequence, 
                bloom, k, hashed_minimizer,
                kmer_number, all_addresses, address_mask, l);
        } else {
            let _ = 
                handle_super_kmer_u128(super_kmers_positions[super_kmers_positions.len()-1], 
                (sequence.len()+1-k as usize) as u32,
                &sequence,
                bloom, k, hashed_minimizer,
                kmer_number, all_addresses, address_mask, l);
        }
    }

    //is here only to prevent optimisations in case no bloom filters
    local_kmer_sum
}

fn handle_super_kmer(start_pos: u32, end_pos: u32, sequence: &PackedSeqVec,
    bloom: &BloomFilter, 
    k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, all_addresses: &mut Vec<usize>, 
    address_mask: usize, l: u16) -> usize {
    let mut last_relevant_index: usize = 0;
    for j in (start_pos as usize)..(end_pos as usize) + (k-l) as usize{
        let smer: PackedSeq = sequence.slice(j..j+l as usize);
        let mut hash: u64 = xorshift_u64(smer.as_u64());


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
    k: u16, hashed_minimizer: u64, 
    mut kmer_number: usize, all_addresses: &mut Vec<usize>, 
    address_mask: usize, l: u16) -> usize {
    let mut last_relevant_index: usize = 0;
    for j in (start_pos as usize)..(end_pos as usize) + (k-l) as usize {
        let smer: PackedSeq = sequence.slice(j..j+l as usize);
        let mut hash: u128 = xorshift_u128(smer.as_u128());


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
    bloom: BloomFilter,
    nb_blocks: usize,
    block_size: usize,
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


    //false negatives and false potitives rates
    print_string += &format!("|{:.3}|{:.3}", false_positive_rate, false_negative_rate);

    //duration of overhead decycling set calculation
    print_string += &format!("|{}", duration_overhead_decycling.as_secs());

    println!("{print_string}");
}

#[cfg(test)]
mod tests {
    use super::{
        handle_sequence, handle_super_kmer, handle_super_kmer_u128, resolve_input_path, BloomFilter, Decycler,
    };
    use packed_seq::{PackedSeqVec, SeqVec};

    fn build_bloom(k: usize) -> BloomFilter {
        BloomFilter::new(1 << 20, 3, k, 1 << 10, 1 << 10)
    }

    fn build_decycler(m: u16) -> Decycler {
        let mut decycler = Decycler::new(m);
        decycler.compute_blocks();
        decycler
    }

    #[test]
    fn handle_sequence_returns_zero_for_short_inputs() {
        let bloom = build_bloom(5);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(b"ACGT");
        let mut addresses = vec![0; 64];

        assert_eq!(
            handle_sequence(&bloom, sequence, 5, 3, 1 << 10, false, &mut addresses, &decycler, 3),
            0
        );
    }

    #[test]
    fn handle_sequence_no_bloom_returns_non_zero_hash_sum() {
        let bloom = build_bloom(5);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let mut addresses = vec![0; 64];

        let sum = handle_sequence(&bloom, sequence, 5, 3, 1 << 10, true, &mut addresses, &decycler, 3);
        assert_ne!(sum, 0);
    }

    #[test]
    fn handle_sequence_populates_bloom_for_u64_queries() {
        let bloom = build_bloom(5);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let mut addresses = vec![0; 64];

        let _ = handle_sequence(&bloom, sequence.clone(), 5, 3, 1 << 10, false, &mut addresses, &decycler, 3);
        let results = bloom.check_sequence(sequence, 5, 3, 3, &decycler);

        assert_eq!(results.len(), 4);
        assert!(results.iter().all(|present| *present));
    }

    #[test]
    fn handle_sequence_populates_bloom_for_u128_queries() {
        let bloom = build_bloom(40);
        let decycler = build_decycler(3);
        let sequence = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );
        let mut addresses = vec![0; 256];

        let _ = handle_sequence(&bloom, sequence.clone(), 40, 3, 1 << 10, false, &mut addresses, &decycler, 32);
        let results = bloom.check_sequence_u128(sequence, 40, 3, 32, &decycler);

        assert_eq!(results.len(), 9);
        assert!(results.iter().all(|present| *present));
    }

    #[test]
    fn handle_super_kmer_returns_updated_counter() {
        let bloom = build_bloom(5);
        let sequence = PackedSeqVec::from_ascii(b"ACGTACGT");
        let mut addresses = vec![0; 64];

        let count = handle_super_kmer(0, 2, &sequence, &bloom, 5, 0, 7, &mut addresses, (1 << 10) - 1, 3);
        assert_eq!(count, 11);
    }

    #[test]
    fn handle_super_kmer_u128_returns_updated_counter() {
        let bloom = build_bloom(40);
        let sequence = PackedSeqVec::from_ascii(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        );
        let mut addresses = vec![0; 256];

        let count = handle_super_kmer_u128(0, 2, &sequence, &bloom, 40, 0, 5, &mut addresses, (1 << 10) - 1, 32);
        assert_eq!(count, 15);
    }

    #[test]
    fn resolve_input_path_returns_positional_input() {
        let path = resolve_input_path(Some("reads.fa"), None);
        assert_eq!(path, "reads.fa");
    }

    #[test]
    fn resolve_input_path_returns_indexed_option_path() {
        let path = resolve_input_path(None, Some("indexed.fa"));
        assert_eq!(path, "indexed.fa");
    }

    #[test]
    #[should_panic(expected = "an input file path is required")]
    fn resolve_input_path_requires_one_input_source() {
        let _ = resolve_input_path(None, None);
    }

    #[test]
    #[should_panic(expected = "use either the positional input or --indexed-file, not both")]
    fn resolve_input_path_rejects_both_input_sources() {
        let _ = resolve_input_path(Some("reads.fa"), Some("indexed.fa"));
    }
}
