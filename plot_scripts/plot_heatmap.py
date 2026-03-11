#!/usr/bin/env python3
"""Create a chr12 cenhap typing accuracy heatmap."""

import argparse # Command line argument parsing
import os # List files in log directory
from typing import Dict, List, Tuple # Type hints

import matplotlib.pyplot as plt # Basic plotting
import numpy as np # Help to make heatmap

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot cenhap typing heatmap.")
    parser.add_argument("--log-dir", required=True)
    parser.add_argument("--cenhap-table", required=True)
    parser.add_argument("--output-file", required=True)

    return parser.parse_args()

def load_truth_table(path: str) -> Dict[str, str]:
    """Read truth cenhap table into {haplo: cenhap} dict."""
    truth = {}
    with open(path) as file:
        file.readline() # header
        for line in file:
            parts = line.strip().split('\t')
            truth[parts[0]] = parts[1]
    return truth

def parse_guess_from_log(path: str) -> str:
    """Find "Best guess" line and return guessed cenhap."""
    with open(path) as file:
        for line in file:
            if line.startswith("Best guess"):
                return line.strip().split()[-1]
    return None

def collect_guesses(log_dir: str) -> Dict[str, str]:
    """Parse all logs to get guesses as {haplo: cenhap} dict."""
    guesses = {}
    for fname in os.listdir(log_dir):
        if not fname.endswith(".log"):
            continue
        hap = fname[:-4]
        path = os.path.join(log_dir, fname)
        guess = parse_guess_from_log(path)
        if guess:
            guesses[hap] = guess
    return guesses

def build_matrix(truth_map: Dict[str, str], guess_map: Dict[str, str]
                 ) -> Tuple[np.ndarray, List[str], List[str]]:
    """Build confusion matrix counts."""
    rows = []
    cols = []

    # Find all truth/guess pairs
    pairs = []
    for hap, truth in truth_map.items():
        guess = guess_map.get(hap)
        if guess is None:
            continue
        pairs.append((truth, guess))
        rows.append(truth)
        cols.append(guess)

    # Build rows/cols
    truth_labels = sorted(set(rows))
    guess_labels = sorted(set(cols))

    row_index = {x: i for i, x in enumerate(truth_labels)}
    col_index = {x: i for i, x in enumerate(guess_labels)}

    # Fill in matrix
    mat = np.zeros((len(truth_labels), len(guess_labels)), dtype=int)
    for truth, guess in pairs:
        mat[row_index[truth], col_index[guess]] += 1

    return mat, truth_labels, guess_labels

def plot_heatmap(counts: np.ndarray, truth_labels: List[str], 
                 guess_labels: List[str], output_file: str) -> None:
    """Plot row-normalized heatmap with count annotations."""
    row_sums = counts.sum(axis=1, keepdims=True)
    norm = np.divide(counts, row_sums, where=row_sums != 0)

    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.imshow(norm, cmap="viridis", vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(guess_labels)))
    ax.set_yticks(np.arange(len(truth_labels)))

    ax.set_xticklabels(guess_labels, rotation=45, ha="right")
    ax.set_yticklabels(truth_labels)

    ax.set_xlabel("Guessed cenhap")
    ax.set_ylabel("Truth cenhap")
    ax.set_title("chr12 cenhap typing")

    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            val = counts[i, j]
            if val > 0:
                ax.text(j, i, str(val), ha="center", va="center", color="white")

    fig.colorbar(im, ax=ax, label="Row-normalized fraction")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

if __name__ == "__main__":
    args = parse_args()
    truth_map = load_truth_table(args.cenhap_table)
    guess_map = collect_guesses(args.log_dir)

    counts, truth_labels, guess_labels = build_matrix(truth_map, guess_map)

    if counts.size == 0:
        raise RuntimeError("No overlapping truth/guess assignments found.")

    plot_heatmap(counts, truth_labels, guess_labels, args.output_file)