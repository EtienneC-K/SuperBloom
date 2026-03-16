# Benchmark Scripts

This folder provides Python scripts to benchmark the `benchmark` Rust binary by sweeping one parameter at a time while keeping all others at default values.

Each sweep script:

- takes an index FASTA/FASTQ and a query FASTA/FASTQ,
- runs multiple benchmark points,
- writes a TSV file,
- writes a PNG plot.

All runs force `MinimizerMode::Simd` via the Rust benchmark binary.

## Prerequisites

- Python 3
- `matplotlib` (`pip install matplotlib`)
- Rust toolchain (`cargo`)

## Available scripts

- `sweep_k.py`
- `sweep_m.py`
- `sweep_s.py`
- `sweep_n_hashes.py`
- `sweep_size_exponent.py`
- `sweep_block_size_exponent.py`
- `sweep_threads.py`
- `run_all_sweeps.py`

## Example: one sweep

```bash
python benchmark/sweep_threads.py \
  --index-fasta data/ecoli.fa.zst \
  --query-fasta data/ecoli.fa.zst
```

Outputs:

- `benchmark/results/sweep_threads.tsv`
- `benchmark/results/sweep_threads.png`

## Example: run all sweeps

```bash
python benchmark/run_all_sweeps.py \
  --index-fasta data/ecoli.fa.zst \
  --query-fasta data/ecoli.fa.zst
```
