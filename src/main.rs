use bloomybloom::{MinimizerMode, SuperBloom, SuperBloomConfig};
use needletail::parse_fastx_file;
use rayon::ThreadPoolBuilder;
use std::error::Error;
use std::fs;
use std::path::PathBuf;

fn first_query_of_len(path: &str, query_len: usize) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut reader = parse_fastx_file(path)?;
    while let Some(record_result) = reader.next() {
        let record = record_result?;
        let seq = record.seq().as_ref().to_vec();
        if seq.len() >= query_len {
            return Ok(seq[..query_len].to_vec());
        }
    }
    Err(format!("no sequence of length >= {query_len} found in {path}").into())
}

fn main() -> Result<(), Box<dyn Error>> {
    println!("SuperBloom Library Showcase");
    println!("===========================");
    println!("Dataset: data/ecoli.fa.zst (index + query source)");

    let available_threads = std::thread::available_parallelism()
        .map(|threads| threads.get())
        .unwrap_or(8);
    let rayon_threads = available_threads.max(8);
    let _ = ThreadPoolBuilder::new()
        .num_threads(rayon_threads)
        .build_global();
    println!("Rayon thread pool configured with {rayon_threads} threads.");

 
    let config = SuperBloomConfig {
        k: 31,
        m: 21,
        s: 27,
        n_hashes: 4,
        size_exponent: 35,
        block_size_exponent: 13,
        minimizer_mode: MinimizerMode::Simd,
    };
    let index_fasta = "data/ecoli.fa.zst";
    let query_fasta = "data/ecoli.fa.zst";
    let query_len = 100usize;
    let query_sequence = first_query_of_len(query_fasta, query_len)?;

    println!("\n1) SuperBloom::new(config)");
    // Create a new mutable index from fully manual parameters.
    let mut bloom = SuperBloom::new(config)?;
    println!("   created index with config: {:?}", bloom.config());

    println!("\n2) SuperBloom::add_sequence(&[u8])");
    // Add one in-memory DNA sequence.
    let added_from_query = bloom.add_sequence(&query_sequence)?;
    println!("   added k-mers from one 100bp query sequence: {added_from_query}");
    println!("   total inserted so far: {}", bloom.inserted_kmers());

    println!("\n3) SuperBloom::add_fasta(path)");
    // Add every record from compressed data/ecoli.fa.zst to the index.
    let add_report = bloom.add_fasta(index_fasta)?;
    println!("   records processed: {}", add_report.records_processed);
    println!("   records indexed:   {}", add_report.records_indexed);
    println!("   k-mers added:      {}", add_report.kmers_added);
    println!("   total inserted now: {}", bloom.inserted_kmers());

    println!("\n4) SuperBloom::query_sequence(&[u8])");
    // Querying auto-switches internally to lock-free query mode.
    let query_hits = bloom.query_sequence(&query_sequence)?;
    let positives = query_hits.iter().filter(|&&hit| hit).count();
    println!("   query length (bp): {}", query_len);
    println!("   windows queried:   {}", query_hits.len());
    println!("   positive windows: {positives}");

    println!("\n5) SuperBloom::query_fasta(path)");
    // Query a whole FASTA file and report aggregate hit counts.
    let query_report = bloom.query_fasta(query_fasta)?;
    println!("   records processed: {}", query_report.records_processed);
    println!("   k-mers queried:    {}", query_report.queried_kmers);
    println!("   positive k-mers:   {}", query_report.positive_kmers);
    println!(
        "   found by filter:   {} / {} k-mers",
        query_report.positive_kmers, query_report.queried_kmers
    );

    println!("\n6) SuperBloom::save(path) + SuperBloom::load(path)");
    // Save current index to disk, then load it back.
    let save_path: PathBuf = std::env::temp_dir().join("superbloom_ecoli_demo.sbf");
    println!("   serializing full  index to disk");
    bloom.save(&save_path)?;
    println!("   serialized frozen index to: {}", save_path.display());

    let mut loaded = SuperBloom::load(&save_path)?;
    println!("   loaded index from file.");
    println!("   loaded inserted_kmers(): {}", loaded.inserted_kmers());

    println!("\n7) Insert after loading");
    // Insertion after loading auto-switches internally back to mutable mode.
    let added_after_load = loaded.add_sequence(&query_sequence)?;
    println!("   added k-mers after load: {added_after_load}");
    let hits_after_insert = loaded.query_sequence(&query_sequence)?;
    let positives_after_insert = hits_after_insert.iter().filter(|&&hit| hit).count();
    println!(
        "   query after insert: {} / {} positive windows",
        positives_after_insert,
        hits_after_insert.len()
    );

    println!("\n8) Re-query after loading");
    // Run the same whole-file query again on the loaded index.
    let loaded_hits = loaded.query_sequence(&query_sequence)?;
    let loaded_positive = loaded_hits.iter().filter(|&&hit| hit).count();
    let loaded_report = loaded.query_fasta(query_fasta)?;
    println!(
        "   query_sequence positives (loaded): {} / {}",
        loaded_positive,
        loaded_hits.len()
    );
    println!(
        "   query_fasta found (loaded): {} / {} k-mers",
        loaded_report.positive_kmers, loaded_report.queried_kmers
    );

    let _ = fs::remove_file(&save_path);

    println!("\nShowcase complete.");
    Ok(())
}
