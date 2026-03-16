use bloomybloom::{FrozenSuperBloom, MinimizerMode, SuperBloom, SuperBloomConfig};
use needletail::parse_fastx_file;
use rayon::ThreadPoolBuilder;
use std::error::Error;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

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
    let total_start = Instant::now();
    println!("SuperBloom Library Showcase");
    println!("===========================");
    println!("Dataset: ecoli.fa (index + query source)");

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
    let index_fasta = "ecoli.fa";
    let query_fasta = "ecoli.fa";
    let query_len = 100usize;
    let pre_query_start = Instant::now();
    let query_sequence = first_query_of_len(query_fasta, query_len)?;
    let pre_query_seconds = pre_query_start.elapsed().as_secs_f64();

    println!("\n1) SuperBloom::new(config)");
    let phase_1_start = Instant::now();
    let mut bloom = SuperBloom::new(config)?;
    println!("   created index with config: {:?}", bloom.config());
    let phase_1_seconds = phase_1_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_1_seconds:.3}s");

    println!("\n2) SuperBloom::add_sequence(&[u8])");
    let phase_2_start = Instant::now();
    let added_from_query = bloom.add_sequence(&query_sequence)?;
    println!("   added k-mers from one 100bp query sequence: {added_from_query}");
    println!("   total inserted so far: {}", bloom.inserted_kmers());
    let phase_2_seconds = phase_2_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_2_seconds:.3}s");

    println!("\n3) SuperBloom::add_fasta(path)");
    let phase_3_start = Instant::now();
    let add_report = bloom.add_fasta(index_fasta)?;
    println!("   records processed: {}", add_report.records_processed);
    println!("   records indexed:   {}", add_report.records_indexed);
    println!("   k-mers added:      {}", add_report.kmers_added);
    println!("   total inserted now: {}", bloom.inserted_kmers());
    let phase_3_seconds = phase_3_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_3_seconds:.3}s");

    println!("\n4) SuperBloom::into_frozen()");
    let phase_4_start = Instant::now();
    let frozen = bloom.into_frozen();
    println!("   frozen inserted_kmers(): {}", frozen.inserted_kmers());
    println!("   frozen config(): {:?}", frozen.config());
    let phase_4_seconds = phase_4_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_4_seconds:.3}s");

    println!("\n5) FrozenSuperBloom::query_sequence(&[u8])");
    let phase_5_start = Instant::now();
    let query_hits = frozen.query_sequence(&query_sequence);
    let positives = query_hits.iter().filter(|&&hit| hit).count();
    println!("   query length (bp): {}", query_len);
    println!("   windows queried:   {}", query_hits.len());
    println!("   positive windows: {positives}");
    let phase_5_seconds = phase_5_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_5_seconds:.3}s");

    println!("\n6) FrozenSuperBloom::query_fasta(path)");
    let phase_6_start = Instant::now();
    let query_report = frozen.query_fasta(query_fasta)?;
    println!("   records processed: {}", query_report.records_processed);
    println!("   k-mers queried:    {}", query_report.queried_kmers);
    println!("   positive k-mers:   {}", query_report.positive_kmers);
    println!(
        "   found by filter:   {} / {} k-mers",
        query_report.positive_kmers, query_report.queried_kmers
    );
    let phase_6_seconds = phase_6_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_6_seconds:.3}s");

    println!("\n7) FrozenSuperBloom::save(path) + FrozenSuperBloom::load(path)");
    let phase_7_start = Instant::now();
    let save_path: PathBuf = std::env::temp_dir().join("superbloom_ecoli_demo.sbf");
    println!("   serializing full 4GB-bit index to disk (this step is I/O heavy)...");
    frozen.save(&save_path)?;
    println!("   serialized frozen index to: {}", save_path.display());

    let loaded = FrozenSuperBloom::load(&save_path)?;
    println!("   loaded index from file.");
    println!("   loaded inserted_kmers(): {}", loaded.inserted_kmers());
    let phase_7_seconds = phase_7_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_7_seconds:.3}s");

    println!("\n8) Re-query after loading");
    let phase_8_start = Instant::now();
    let loaded_hits = loaded.query_sequence(&query_sequence);
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
    let phase_8_seconds = phase_8_start.elapsed().as_secs_f64();
    println!("   phase duration: {phase_8_seconds:.3}s");

    let _ = fs::remove_file(&save_path);

    let total_seconds = total_start.elapsed().as_secs_f64();
    println!("\nPhase timing summary:");
    println!("   query extraction (100bp): {pre_query_seconds:.3}s");
    println!("   1) new(config):            {phase_1_seconds:.3}s");
    println!("   2) add_sequence:           {phase_2_seconds:.3}s");
    println!("   3) add_fasta:              {phase_3_seconds:.3}s");
    println!("   4) into_frozen:            {phase_4_seconds:.3}s");
    println!("   5) query_sequence:         {phase_5_seconds:.3}s");
    println!("   6) query_fasta:            {phase_6_seconds:.3}s");
    println!("   7) save+load:              {phase_7_seconds:.3}s");
    println!("   8) re-query after load:    {phase_8_seconds:.3}s");
    println!("   total runtime:             {total_seconds:.3}s");

    println!("\nShowcase complete.");
    Ok(())
}
