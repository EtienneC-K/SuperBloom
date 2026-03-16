#!/usr/bin/env python3
"""Run all parameter sweeps one after another."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Run all benchmark sweeps")
    parser.add_argument("--index-fasta", required=True)
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--output-dir", default=str(Path(__file__).resolve().parent / "results"))
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    scripts = [
        "sweep_k.py",
        "sweep_m.py",
        "sweep_s.py",
        "sweep_n_hashes.py",
        "sweep_size_exponent.py",
        "sweep_block_size_exponent.py",
        "sweep_threads.py",
    ]

    for script in scripts:
        script_path = script_dir / script
        cmd = [
            sys.executable,
            str(script_path),
            "--index-fasta",
            args.index_fasta,
            "--query-fasta",
            args.query_fasta,
            "--output-dir",
            args.output_dir,
        ]
        print(f"\n==> Running {script}")
        subprocess.run(cmd, check=True)

    print("\nAll sweeps completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
