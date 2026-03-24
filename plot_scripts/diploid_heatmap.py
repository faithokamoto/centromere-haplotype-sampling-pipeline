#!/usr/bin/env python3

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot heatmap of diploid cenhap accuracy."
    )
    parser.add_argument(
        "--input-tsv",
        required=True,
        help="TSV from evaluation script"
    )
    parser.add_argument(
        "--output",
        default="heatmap.png",
        help="Output image file (default: heatmap.png)"
    )
    return parser.parse_args()


def canonicalize_pair(pair_str):
    """
    Convert '1,4' or '4,1' -> '1/4' (sorted, canonical form)
    """
    parts = [p.strip() for p in pair_str.split(",")]
    parts_sorted = sorted(parts, key=lambda x: (int(x) if x.isdigit() else x))
    return "/".join(parts_sorted)


def load_rows(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append({
                "truth": canonicalize_pair(r["true_cenhaps"]),
                "guess": canonicalize_pair(r["guessed_cenhaps"])
            })
    return rows


def make_heatmap(ax, rows):
    # UNION of all labels → ensures square matrix with same axes
    labels = sorted({r['truth'] for r in rows} | {r['guess'] for r in rows})
    index = {v: i for i, v in enumerate(labels)}

    counts = np.zeros((len(labels), len(labels)))

    for r in rows:
        counts[index[r['truth']], index[r['guess']]] += 1

    # Normalize row-wise
    norm = counts.copy()
    for i in range(len(labels)):
        s = counts[i].sum()
        if s > 0:
            norm[i] /= s

    im = ax.imshow(norm, cmap='viridis', vmin=0, vmax=1)

    # annotate counts
    for i in range(len(labels)):
        for j in range(len(labels)):
            if counts[i, j] > 0:
                ax.text(
                    j, i, int(counts[i, j]),
                    ha='center', va='center',
                    backgroundcolor='white', color='black'
                )

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)

    ax.set_xlabel('Guessed cenhap pair')
    ax.set_ylabel('True cenhap pair')

    ax.set_box_aspect(1)

    return im


if __name__ == "__main__":
    args = parse_args()

    rows = load_rows(args.input_tsv)

    fig, ax = plt.subplots(figsize=(8, 8))
    make_heatmap(ax, rows)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()
