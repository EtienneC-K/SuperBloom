# SuperBloom

![tests](https://github.com/EtienneC-K/SuperBloom/workflows/tests/badge.svg)

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
cargo run -r --bin superbloom
```

### Run the benchmark binary

```bash
cargo run --release --bin benchmark -- \
  --index-fasta data/ecoli.fa.zst \
  --query-fasta data/ecoli.fa.zst \
  --k 31 --m 21 --s 27 \
  --n-hashes 8 \
  --size-exponent 35 \
  --block-size-exponent 9 \
  --threads 8
```

`benchmark` always uses `MinimizerMode::Simd` (other minimizer modes are experimental).

For automated parameter sweeps (TSV + plots), see [benchmark/README.md](/home/nadine/Code/SuperBloom/benchmark/README.md).

The showcase:

1. builds an index
2. inserts one 100bp query sequence
3. inserts all records from `data/ecoli.fa.zst`
4. queries sequence + FASTA
5. saves to disk
6. reloads
7. inserts again after loading,
8. re-queries.

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
- `set_threads(usize) -> Result<(), SuperBloomError>`
- `clear_threads()`, `threads()`
- `inserted_kmers()`, `config()`



## Fully Commented Example (Library Usage)

```rust
use superbloom::{MinimizerMode, SuperBloom, SuperBloomConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1) Start from crate defaults.
    // Defaults are: k=31, m=21, s=27, h=8, bit_vector_size_exponent=35 (4 GiB),
    // block_size_exponent=9 (512 bits), minimizer_mode=Simd.
    let config = SuperBloomConfig::default();

    // 2) Build mutable index.
    let mut sb = SuperBloom::new(config)?;

    // 3) Configure threads used by parallel FASTA indexing.
    sb.set_threads(8)?;

    // 4) Insert one in-memory DNA sequence.
    let added = sb.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;
    println!("added {added} k-mers from memory");

    // 5) Insert compressed FASTA/FASTQ from disk.
    // needletail auto-detects compression (gz/bz2/xz/zst) and format (FASTA/FASTQ).
    let add_report = sb.add_fasta("data/ecoli.fa.zst")?;
    println!("indexed {} records", add_report.records_indexed);

    // 6) Query from memory (one bool per k-mer window).
    let hits = sb.query_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;
    let positives = hits.iter().filter(|&&x| x).count();
    println!("memory query positives: {positives}/{}", hits.len());

    // 7) Query an entire FASTA/FASTQ.
    let q = sb.query_fasta("data/ecoli.fa.zst")?;
    println!("fasta positives: {}/{}", q.positive_kmers, q.queried_kmers);

    // 8) Save and load.
    sb.save("/tmp/demo.sbf")?;
    let mut sb2 = SuperBloom::load("/tmp/demo.sbf")?;

    // 9) Insert after querying is allowed.
    // Internally, the index auto-thaws/freeze as needed.
    let _ = sb2.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;

    Ok(())
}
```

## Parameter Suggestions

All geometry is explicit in `SuperBloomConfig`, and threading is controlled at runtime on `SuperBloom`.

### Defaults

`SuperBloomConfig::default()` uses:

- `k = 31`
- `m = 21`
- `s = 27` (`k - 4`)
- `n_hashes = 8`
- `bit_vector_size_exponent = 35` (4 GiB bit-array)
- `block_size_exponent = 9` (512-bit blocks)
- `minimizer_mode = MinimizerMode::Simd`

Runtime defaults:

- `SuperBloom::new(...)` starts with `8` indexing threads for `add_fasta`.
- Override with `set_threads(n)`, reset with `clear_threads()`.

### Meaning and Influence

- `k` (k-mer length, default `31`)
  - Higher `k`: usually more specific matching and fewer random hits, but more sensitivity to sequencing errors/variants.
  - For fixed `m`, higher `k` also tends to create longer super-k-mers, which improves streaming locality and can improve performance.
  - Lower `k`: more tolerant matching, but higher chance of ambiguous/false matches.

- `m` (minimizer length, default `21`, must satisfy `m < k`)
  - Higher `m`: generally more minimizer changes and less super-k-mer grouping.
  - Lower `m`: larger super-k-mer groups and stronger locality, but potentially more collisions/less specificity.

- `s` (findere-like subword length, default `27`)
  - Rule of thumb: `s = k - 4` is a strong baseline.
  - Higher `s` (closer to `k`): behavior gets closer to plain k-mer checking.
  - Lower `s`: stronger false-positive suppression from overlap consistency, but more subword checks per k-mer.

- `n_hashes` (default `8`)
  - Higher values can reduce false positives up to a point.
  - Too high harms speed (more hash probes and memory touches).
  - Too low is faster but generally increases false positives.

- `bit_vector_size_exponent` (total filter bits = `2^bit_vector_size_exponent`, default `35`)
  - Main memory/accuracy knob.
  - Higher values use more RAM and usually lower false positives.
  - Lower values save RAM but raise false positives under the same workload.
  - RAM bytes used by bit-array: `2^(bit_vector_size_exponent - 3)`.
  - Examples: `30` = 128 MiB, `33` = 1 GiB, `35` = 4 GiB.

- `block_size_exponent` (block bits = `2^block_size_exponent`, default `9`)
  - Controls locality granularity.
  - Smaller blocks improve cache fit/locality but can increase block pressure/collisions.
  - Larger blocks reduce per-block pressure but may reduce cache efficiency.

- `minimizer_mode` (default `Simd`)
  - `Simd`: recommended stable mode.
  - `Decycling` and `OpenClosed { t }`: experimental; library emits a warning when selected.

- `threads` (runtime, default `8`)
  - Set with `set_threads(n)`.
  - Controls parallelism for `add_fasta` and `query_fasta`.
  - More threads usually speed up those operations on larger inputs, up to memory bandwidth / CPU limits.

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
- Large `bit_vector_size_exponent` values (e.g., `35`) require significant RAM.
- Large number of hash function or large block size will harm perfomances





## Citation

To be done
