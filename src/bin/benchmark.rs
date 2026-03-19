use clap::{Arg, Command, value_parser};
use std::error::Error;
use std::time::Instant;
use superbloom::{MinimizerMode, SuperBloom, SuperBloomConfig};

// TODO use Cli struct
fn cli() -> Command {
    Command::new("superbloom-benchmark")
        .about("Benchmark-oriented SuperBloom runner (SIMD minimizers only)")
        .arg(
            Arg::new("index-fasta")
                .long("index-fasta")
                .short('i')
                .value_name("PATH")
                .required(true)
                .help("Input FASTA/FASTQ file to index (compressed formats supported)"),
        )
        .arg(
            Arg::new("query-fasta")
                .long("query-fasta")
                .short('q')
                .value_name("PATH")
                .required(true)
                .help("Input FASTA/FASTQ file to query (compressed formats supported)"),
        )
        .arg(
            Arg::new("k")
                .long("k")
                .value_parser(value_parser!(u16))
                .default_value("31")
                .help("k-mer length"),
        )
        .arg(
            Arg::new("m")
                .long("m")
                .value_parser(value_parser!(u16))
                .default_value("21")
                .help("minimizer length (must be <= k and < 32)"),
        )
        .arg(
            Arg::new("s")
                .long("s")
                .value_parser(value_parser!(u16))
                .help("s-mer length used by findere-like checks (default: k-4)"),
        )
        .arg(
            Arg::new("n-hashes")
                .long("n-hashes")
                .value_parser(value_parser!(usize))
                .default_value("8")
                .help("number of hash probes"),
        )
        .arg(
            Arg::new("size-exponent")
                .long("size-exponent")
                .value_parser(value_parser!(u8))
                .default_value("35")
                .help("bit-array size exponent: total bits = 2^size_exponent"),
        )
        .arg(
            Arg::new("block-size-exponent")
                .long("block-size-exponent")
                .value_parser(value_parser!(u8))
                .default_value("9")
                .help("block size exponent: block bits = 2^block_size_exponent"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .value_parser(value_parser!(usize))
                .default_value("8")
                .help("thread count used for add_fasta/query_fasta"),
        )
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = cli().get_matches();

    let index_fasta = matches
        .get_one::<String>("index-fasta")
        .expect("required argument index-fasta")
        .as_str();
    let query_fasta = matches
        .get_one::<String>("query-fasta")
        .expect("required argument query-fasta")
        .as_str();
    let k = *matches.get_one::<u16>("k").expect("defaulted");
    let m = *matches.get_one::<u16>("m").expect("defaulted");
    let s = matches
        .get_one::<u16>("s")
        .copied()
        .unwrap_or_else(|| k.saturating_sub(4).max(1));
    let n_hashes = *matches.get_one::<usize>("n-hashes").expect("defaulted");
    let size_exponent = *matches.get_one::<u8>("size-exponent").expect("defaulted");
    let block_size_exponent = *matches
        .get_one::<u8>("block-size-exponent")
        .expect("defaulted");
    let threads = *matches.get_one::<usize>("threads").expect("defaulted");

    let config = SuperBloomConfig {
        k,
        m,
        s,
        n_hashes,
        size_exponent,
        block_size_exponent,
        minimizer_mode: MinimizerMode::Simd, // non-SIMD modes are experimental
    };

    println!("SuperBloom Benchmark");
    println!("====================");
    println!("index fasta: {index_fasta}");
    println!("query fasta: {query_fasta}");
    println!("config: {:?}", config);
    println!("threads: {threads}");
    println!("minimizer_mode forced to SIMD");

    let total_start = Instant::now();

    let create_start = Instant::now();
    let mut bloom = SuperBloom::new(config)?;
    bloom.set_threads(threads)?;
    let create_secs = create_start.elapsed().as_secs_f64();
    println!("\n1) build index object: {create_secs:.3}s");

    let index_start = Instant::now();
    let add_report = bloom.add_fasta(index_fasta)?;
    let index_secs = index_start.elapsed().as_secs_f64();
    println!("2) index fasta: {index_secs:.3}s");
    println!("   records processed: {}", add_report.records_processed);
    println!("   records indexed:   {}", add_report.records_indexed);
    println!("   k-mers added:      {}", add_report.kmers_added);
    if index_secs > 0.0 {
        println!(
            "   indexing throughput: {:.0} k-mers/s",
            add_report.kmers_added as f64 / index_secs
        );
    }

    let query_start = Instant::now();
    let query_report = bloom.query_fasta(query_fasta)?;
    let query_secs = query_start.elapsed().as_secs_f64();
    println!("3) query fasta: {query_secs:.3}s");
    println!("   records processed: {}", query_report.records_processed);
    println!("   queried k-mers:    {}", query_report.queried_kmers);
    println!("   positive k-mers:   {}", query_report.positive_kmers);
    if query_secs > 0.0 {
        println!(
            "   query throughput:   {:.0} k-mers/s",
            query_report.queried_kmers as f64 / query_secs
        );
    }

    let total_secs = total_start.elapsed().as_secs_f64();
    println!("\nTotal: {total_secs:.3}s");
    Ok(())
}
