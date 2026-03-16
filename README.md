# SuperBloom

A Rust implementation of a Super Bloom filter for streaming DNA `k`-mer indexing and querying.


## What This Crate Provides

- A **library API** to:
  - build a SuperBloom index from explicit parameters,
  - insert DNA from memory (`add_sequence`) or from FASTA/FASTQ (`add_fasta`),
  - query from memory (`query_sequence`) or FASTA/FASTQ (`query_fasta`),
  - serialize/deserialize indexes (`save` / `load`).
- A **single showcase binary** (`src/main.rs`) that demonstrates the full API on the bundled compressed E. coli reference genome.

## Repository Layout

- `src/lib.rs`: public library API (`SuperBloom`, `SuperBloomConfig`, reports, errors)
- `src/bloom.rs`: low-level Bloom/SuperBloom data structures and query/insert kernels
- `src/main.rs`: end-to-end showcase executable
- `data/ecoli.fa.zst`: compressed example genome used by showcase
- `data/Superbloom.pdf`: paper used as conceptual reference

## Build, Test, Run

### Prerequisites

- Rust toolchain (stable)
- Standard build tools for crates in `Cargo.toml`

### Build

```bash
cargo build -r
```

### Run tests

```bash
cargo test
```

### Run the showcase binary

```bash
cargo run -r
```

The showcase:

1. builds an index
2. inserts one 100bp query sequence
3. inserts all records from `data/ecoli.fa.zst`
4. queries sequence + FASTA
5. saves to disk
6. reload
7. re-queries
8. inserts again after queries

## Library API 

Main public types:

- `SuperBloomConfig`
- `SuperBloom`
- `FrozenSuperBloom`
- `AddReport`
- `QueryReport`
- `SuperBloomError`
- `MinimizerMode`

### Core `SuperBloom` methods

- `SuperBloom::new(config)`
- `add_sequence(&[u8]) -> Result<u64, SuperBloomError>`
- `add_fasta(path) -> Result<AddReport, SuperBloomError>`
- `query_sequence(&[u8]) -> Result<Vec<bool>, SuperBloomError>`
- `query_fasta(path) -> Result<QueryReport, SuperBloomError>`
- `save(path) -> Result<(), SuperBloomError>`
- `load(path) -> Result<SuperBloom, SuperBloomError>`
- `inserted_kmers()`, `config()`



## Fully Commented Example (Library Usage)

```rust
use bloomybloom::{MinimizerMode, SuperBloom, SuperBloomConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1) Configure geometry and algorithmic parameters explicitly.
    let config = SuperBloomConfig {
        k: 31,
        m: 21,
        s: 27,                    // findere-like setting: s = k - 4
        n_hashes: 4,
        size_exponent: 33,        // 2^33 bits = 1 GiB of RAM for the bit-array
        block_size_exponent: 13,  // 8192 bits per block
        minimizer_mode: MinimizerMode::Simd,
    };

    // 2) Build mutable index.
    let mut sb = SuperBloom::new(config)?;

    // 3) Insert one in-memory DNA sequence.
    let added = sb.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;
    println!("added {added} k-mers from memory");

    // 4) Insert compressed FASTA/FASTQ from disk.
    // needletail auto-detects compression (gz/bz2/xz/zst) and format (FASTA/FASTQ).
    let add_report = sb.add_fasta("data/ecoli.fa.zst")?;
    println!("indexed {} records", add_report.records_indexed);

    // 5) Query from memory (one bool per k-mer window).
    let hits = sb.query_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;
    let positives = hits.iter().filter(|&&x| x).count();
    println!("memory query positives: {positives}/{}", hits.len());

    // 6) Query an entire FASTA/FASTQ.
    let q = sb.query_fasta("data/ecoli.fa.zst")?;
    println!("fasta positives: {}/{}", q.positive_kmers, q.queried_kmers);

    // 7) Save and load.
    sb.save("/tmp/demo.sbf")?;
    let mut sb2 = SuperBloom::load("/tmp/demo.sbf")?;

    // 8) Insert after querying is allowed.
    // Internally, the index auto-thaws/freeze as needed.
    let _ = sb2.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;

    Ok(())
}
```

## Parameter Suggestions

All geometry is explicit in `SuperBloomConfig`.

- `k`: k-mer length.
  - Typical genomic default: `31`.
- `m`: minimizer length (`m < k`).
  - Common fast/accurate setting here: `m=21`.
- `s`: subword length used by findere-like logic (`s <= k`).
  - Good default: `s = k - 4`.
  - For `k=31`, use `s=27`.
- `n_hashes`: number of hash probes per subword.
  - Typical: `3..8`.
  - Good default: `4`.
- `size_exponent`: total bit-array size as `2^size_exponent` bits.
  - RAM bytes for bit-array = `2^(size_exponent - 3)`.
  - Examples:
    - `30` -> 128 MiB
    - `33` -> 1 GiB
    - `35` -> 4 GiB
- `block_size_exponent`: block size as `2^block_size_exponent` bits.
  - Practical range: `12..14`.
  - Current showcase: `13` (8192 bits).
- `minimizer_mode`:
  - `Simd` for speed-focused default,
  - `Decycling` for decycler-based selection,
  - `OpenClosed { t }` for open/closed minimizer strategy.

## What Happens Under the Hood


1. **Blocked layout for locality**  
   The filter is divided into blocks so accesses for related queries stay localized.

2. **Minimizer-based super-k-mer grouping**  
   Consecutive k-mers often share a minimizer. SuperBloom maps these grouped k-mers to the same block, amortizing random accesses over sequence streaming.

3. **findere-style consistency through `s`-mers**  
   Instead of relying only on raw k-mer membership, overlapping subword evidence (`s`-mers) reduces false positives while keeping throughput high.


4. **Parallel ingestion**  
   FASTA insertion is batched and parallelized with Rayon (`PAR_BATCH_RECORDS`)

5. **Compressed FASTA/FASTQ support**  
   Input readers use `needletail::parse_fastx_file`, which supports compressed files (`.gz`, `.bz2`, `.xz`, `.zst`) and auto-detects FASTA/FASTQ.

## Serialization

- Binary format with magic header: `SBLOOM01`.
- `save` serializes config + inserted count + frozen filter shards.
- `load` restores the full structure and keeps it query-ready.


## Performance Notes

- Querying full sequences benefits most (streaming locality, super-k-mer reuse).
- Isolated random k-mer queries behave closer to blocked Bloom behavior.
- Large `size_exponent` values (e.g., `35`) require significant RAM.
- Large number of hash function or large block size will harm perfomances





## Citation

To be done
