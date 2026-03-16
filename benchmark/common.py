#!/usr/bin/env python3
"""Shared utilities for SuperBloom benchmark parameter sweeps."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path
from typing import Any, Callable

DEFAULT_PARAMS: dict[str, int] = {
    "k": 31,
    "m": 21,
    "s": 27,
    "n_hashes": 8,
    "size_exponent": 35,
    "block_size_exponent": 9,
    "threads": 8,
}

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "results"

METRIC_PATTERNS: dict[str, str] = {
    "build_s": r"1\)\s+build index object:\s+([0-9.]+)s",
    "index_s": r"2\)\s+index fasta:\s+([0-9.]+)s",
    "query_s": r"3\)\s+query fasta:\s+([0-9.]+)s",
    "total_s": r"Total:\s+([0-9.]+)s",
    "index_kmers_per_s": r"indexing throughput:\s+([0-9.]+)\s+k-mers/s",
    "query_kmers_per_s": r"query throughput:\s+([0-9.]+)\s+k-mers/s",
}


def parse_values(values_csv: str, caster: Callable[[str], Any]) -> list[Any]:
    values: list[Any] = []
    for raw in values_csv.split(","):
        item = raw.strip()
        if not item:
            continue
        values.append(caster(item))
    if not values:
        raise ValueError("values list is empty")
    return values


def run_benchmark(index_fasta: str, query_fasta: str, params: dict[str, int]) -> str:
    cmd = [
        "cargo",
        "run",
        "--quiet",
        "--release",
        "--bin",
        "benchmark",
        "--",
        "--index-fasta",
        index_fasta,
        "--query-fasta",
        query_fasta,
        "--k",
        str(params["k"]),
        "--m",
        str(params["m"]),
        "--s",
        str(params["s"]),
        "--n-hashes",
        str(params["n_hashes"]),
        "--size-exponent",
        str(params["size_exponent"]),
        "--block-size-exponent",
        str(params["block_size_exponent"]),
        "--threads",
        str(params["threads"]),
    ]
    proc = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "benchmark command failed\n"
            f"command: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )
    return proc.stdout


def parse_metrics(stdout: str) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for key, pattern in METRIC_PATTERNS.items():
        match = re.search(pattern, stdout)
        if not match:
            raise RuntimeError(
                f"failed to parse metric '{key}' from benchmark output:\n{stdout}"
            )
        metrics[key] = float(match.group(1))
    return metrics


def ensure_plotting() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required to generate plots. Install with: pip install matplotlib"
        ) from exc
    return plt


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_results(path: Path, parameter: str, rows: list[dict[str, Any]]) -> None:
    plt = ensure_plotting()

    xs = [row["value"] for row in rows]
    build = [row["build_s"] for row in rows]
    index = [row["index_s"] for row in rows]
    query = [row["query_s"] for row in rows]
    total = [row["total_s"] for row in rows]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(xs, build, marker="o", label="build_s")
    ax.plot(xs, index, marker="o", label="index_s")
    ax.plot(xs, query, marker="o", label="query_s")
    ax.plot(xs, total, marker="o", label="total_s")
    ax.set_xlabel(parameter)
    ax.set_ylabel("seconds")
    ax.set_title(f"SuperBloom sweep: {parameter}")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def run_parameter_sweep_cli(
    *,
    parameter: str,
    default_values: list[int],
    description: str,
    value_caster: Callable[[str], Any] = int,
    link_s_to_k: bool = False,
) -> None:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--index-fasta", required=True, help="Input FASTA/FASTQ to index")
    parser.add_argument("--query-fasta", required=True, help="Input FASTA/FASTQ to query")
    parser.add_argument(
        "--values",
        default=",".join(str(x) for x in default_values),
        help=(
            "Comma-separated values to sweep over for this parameter "
            f"(default: {','.join(str(x) for x in default_values)})"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Output directory for TSV/PNG (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--print-command-output",
        action="store_true",
        help="Print raw benchmark binary output for each run",
    )
    args = parser.parse_args()

    values = parse_values(args.values, value_caster)
    output_dir = Path(args.output_dir)

    rows: list[dict[str, Any]] = []
    for value in values:
        params = dict(DEFAULT_PARAMS)
        params[parameter] = int(value)
        if link_s_to_k:
            params["s"] = max(1, int(params["k"]) - 4)

        stdout = run_benchmark(args.index_fasta, args.query_fasta, params)
        metrics = parse_metrics(stdout)

        if args.print_command_output:
            print(stdout)

        row: dict[str, Any] = {
            "parameter": parameter,
            "value": value,
            **params,
            **metrics,
        }
        rows.append(row)
        print(
            f"{parameter}={value}: "
            f"index={metrics['index_s']:.3f}s "
            f"query={metrics['query_s']:.3f}s "
            f"total={metrics['total_s']:.3f}s"
        )

    fieldnames = [
        "parameter",
        "value",
        "k",
        "m",
        "s",
        "n_hashes",
        "size_exponent",
        "block_size_exponent",
        "threads",
        "build_s",
        "index_s",
        "query_s",
        "total_s",
        "index_kmers_per_s",
        "query_kmers_per_s",
    ]

    tsv_path = output_dir / f"sweep_{parameter}.tsv"
    png_path = output_dir / f"sweep_{parameter}.png"

    write_tsv(tsv_path, rows, fieldnames)
    plot_results(png_path, parameter, rows)

    print(f"\nWrote TSV:  {tsv_path}")
    print(f"Wrote plot: {png_path}")

