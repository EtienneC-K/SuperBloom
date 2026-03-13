#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
import subprocess
from datetime import datetime
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit("matplotlib is required: pip install matplotlib") from exc


BENCHMARK_NAME = "bench_s_vs_ffp"
THREAD_VALUE = 16 #TODO: put correct value for actual bench
K_VALUE = 31
S_VALUES = [20, 26, 28, 30, 31]
H_VALUES = [1, 2, 4, 6, 12]
BLOCK_SIZE = 13
RAM_GB = 2 #TODO: put correct value for actual bench
M_VALUE = 13
REPEATS = 1
BUILD_FIRST = True
USE_INDEXED_FILE_FLAG = False
EXTRA_ARGS: list[str] = []

METRIC_PATTERNS = {
    "queried": r"Number of kmer queried : ([0-9eE+.\-]+)",
    "positives": r"Number of positives : ([0-9eE+.\-]+)",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("index_file")
    parser.add_argument("query_file")
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def results_dir() -> Path:
    path = Path(__file__).resolve().parent / "results"
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_release(root: Path) -> None:
    if BUILD_FIRST:
        subprocess.run(["cargo", "build", "-r"], cwd=root, check=True)


def run_bloom(root: Path, command: list[str]) -> dict[str, float]:
    #to check for false positives
    command.append("--counting")
    completed = subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        capture_output=True,
    )
    metrics: dict[str, float] = {}
    for key, pattern in METRIC_PATTERNS.items():
        match = re.search(pattern, completed.stdout)
        if match is None:
            raise RuntimeError(
                f"Missing metric {key} in output.\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
            )
        metrics[key] = float(match.group(1))
    return metrics


def build_command(index_file: str, query_file: str, s: int, h: int) -> list[str]:
    command = [
        "./target/release/bloomybloom",
        "--query-file",
        str(Path(query_file).expanduser()),
        "--ram",
        str(RAM_GB),
        "--threads",
        str(THREAD_VALUE),
        "-k",
        str(K_VALUE),
        "-m",
        str(M_VALUE),
        "--block-size",
        str(BLOCK_SIZE),
        "-s",
        str(s),
        "--n-hashes",
        str(h),
    ]
    if USE_INDEXED_FILE_FLAG:
        command.extend(["--indexed-file", str(Path(index_file).expanduser())])
    else:
        command.append(str(Path(index_file).expanduser()))
    command.extend(EXTRA_ARGS)
    return command


def aggregate(metrics_list: list[dict[str, float]]) -> dict[str, float]:
    return {
        key: statistics.fmean(run[key] for run in metrics_list)
        for key in METRIC_PATTERNS
    }


def write_tsv(rows: list[dict[str, object]], output_path: Path) -> None:
    fieldnames = [
        "benchmark",
        "index_file",
        "query_file",
        "k",
        "ram_gb",
        "threads",
        "repeats",
        "m",
        "ffp",
        "queried",
        "positives",
        "s",
        "n-hashes",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows: list[dict[str, object]], output_path: Path) -> None:
    x_values = [int(row["s"]) for row in rows]
    queried = [float(row["queried"]) for row in rows]
    positives = [float(row["positives"]) for row in rows]
    ffp = [positives[i]/queried[i] for i in range(len(queried))]


    fig, axes = plt.subplots(1, 1, figsize=(12, 5), sharex=True)

    axes.plot(x_values, ffp, marker="o", label="friendly false positives")
    axes.set_title("False positve rates")
    axes.set_xlabel("s")
    axes.set_ylabel("FFP rate")
    axes.grid(True, alpha=0.3)
    axes.legend()

    fig.suptitle("bloomybloom benchmark: s sweep vs ffp")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    root = repo_root()
    build_release(root)

    rows: list[dict[str, object]] = []
    for i in range (len(S_VALUES)):
        s = S_VALUES[i]
        h = H_VALUES[i]
        metrics_list = [
            run_bloom(root, build_command(args.index_file, args.query_file, s, h))
            for _ in range(REPEATS)
        ]
        metrics = aggregate(metrics_list)
        rows.append({
            "benchmark": BENCHMARK_NAME,
            "index_file": str(Path(args.index_file).expanduser()),
            "query_file": str(Path(args.query_file).expanduser()),
            "k": K_VALUE,
            "ram_gb": RAM_GB,
            "threads": THREAD_VALUE,
            "repeats": REPEATS,
            "m": M_VALUE,
            "s": s,
            "n-hashes": h,
            **metrics,
        })

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = results_dir()
    tsv_path = out_dir / f"{BENCHMARK_NAME}-{timestamp}.tsv"
    png_path = out_dir / f"{BENCHMARK_NAME}-{timestamp}.png"
    write_tsv(rows, tsv_path)
    plot_rows(rows, png_path)
    print(tsv_path)
    print(png_path)


if __name__ == "__main__":
    main()
