use std::error::Error;
use std::time::Instant;
use superbloom::{MinimizerMode, SuperBloom, SuperBloomConfig};

use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "superbloom-benchmark",
    about = "Benchmark-oriented SuperBloom runner (SIMD minimizers only)"
)]
struct Cli {
    /// Input FASTA/FASTQ file to index (compressed formats supported)
    #[arg(long, short = 'i', value_name = "PATH")]
    index_fasta: String,

    /// Input FASTA/FASTQ file to query (compressed formats supported)
    #[arg(long, short = 'q', value_name = "PATH")]
    query_fasta: String,

    /// k-mer length
    #[arg(long, default_value_t = 31)]
    k: u16,

    /// minimizer length (must be <= k and < 32)
    #[arg(long, default_value_t = 21)]
    m: u16,

    /// s-mer length used by findere-like checks (default: k-4)
    #[arg(long)]
    s: Option<u16>,

    /// number of hash probes
    #[arg(long, default_value_t = 8)]
    n_hashes: usize,

    /// bit-array size exponent: total bits = 2^size_exponent
    #[arg(long, default_value_t = 35)]
    size_exponent: u8,

    /// block size exponent: block bits = 2^block_size_exponent
    #[arg(long, default_value_t = 9)]
    block_size_exponent: u8,

    /// thread count used for add_fasta/query_fasta
    #[arg(long, short = 't', default_value_t = 8)]
    threads: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let Cli {
        index_fasta,
        query_fasta,
        k,
        m,
        s,
        n_hashes,
        size_exponent,
        block_size_exponent,
        threads,
    } = Cli::parse();

    let s = s.unwrap_or_else(|| k.saturating_sub(4).max(1));

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
