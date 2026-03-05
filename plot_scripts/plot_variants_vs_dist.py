#!/usr/bin/env python3
"""Plot variant counts vs distance to CHM13 across multiple paths.

Modified from original single-sample plotting script.

This version:
- Accepts a directory of logfiles (<dir>/<path>.log)
- Accepts a file listing path names (one per line)
- Accepts a distance matrix CSV (no header):
      path1,path2,distance
  Only rows where path1 == "CHM13.0" are used.
- Extracts variant counts for haplotypes 1–8
- Produces 8 separate 2x2 scatterplots (one per sampled haplotype)
- Adds linear regression lines (least squares fit)
- Output files are named:
      <prefix>.n1.png
      <prefix>.n2.png
      ...
      <prefix>.n8.png

Each path contributes a single point per panel.

Updated by ChatGPT.
"""

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class VariantCounts:
    chm13_snvs: int
    chm13_svs: int
    alt_snvs: int
    alt_svs: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot variant counts vs distance to CHM13."
    )
    parser.add_argument(
        "--log-dir",
        required=True,
        help="Directory containing <path>.log files."
    )
    parser.add_argument(
        "--path-list",
        required=True,
        help="File with path names (one per line)."
    )
    parser.add_argument(
        "--dist-matrix",
        required=True,
        help="CSV (no header): path1,path2,distance"
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for output plots."
    )
    return parser.parse_args()


def load_distances(dist_matrix: str) -> Dict[str, float]:
    """Return {path2: distance_to_CHM13}."""
    dist_dict = {}

    with open(dist_matrix) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) != 3:
                continue
            path1, path2, dist = parts
            if path1 == "CHM13.0":
                dist_dict[path2] = float(dist)

    return dist_dict


def parse_logfile(logfile: str) -> List[VariantCounts]:
    """Extract VariantCounts for haplotypes 1–8."""
    stats = []

    with open(logfile) as f:
        lines = f.readlines()

    for i in range(len(lines)):
        if lines[i].strip().endswith("SVs on CHM13"):
            chm13_snvs = int(lines[i].split()[0])
            chm13_svs = int(lines[i].split()[3])
            alt_snvs = int(lines[i+1].split()[0])
            alt_svs = int(lines[i+1].split()[3])

            stats.append(
                VariantCounts(
                    chm13_snvs,
                    chm13_svs,
                    alt_snvs,
                    alt_svs
                )
            )

    if len(stats) < 8:
        raise ValueError(f"{logfile} does not contain 8 haplotype runs.")

    return stats[:8]


def add_regression(ax, x, y):
    """Add least-squares regression line to axis."""
    if len(x) < 2:
        return

    m, b = np.polyfit(x, y, 1)
    x_sorted = np.sort(x)
    y_fit = m * x_sorted + b
    ax.plot(x_sorted, y_fit)


def plot_for_haplotype(
    n: int,
    distances: List[float],
    counts: List[VariantCounts],
    output_prefix: str
) -> None:

    fig, axs = plt.subplots(2, 2, figsize=(5, 5))

    x = np.array(distances)

    # Top left: CHM13 SNV
    y = np.array([c.chm13_snvs for c in counts])
    axs[0][0].scatter(x, y)
    add_regression(axs[0][0], x, y)
    axs[0][0].set_title("CHM13 SNVs")
    axs[0][0].set_ylabel("# variants")

    # Top right: CHM13 SV
    y = np.array([c.chm13_svs for c in counts])
    axs[0][1].scatter(x, y)
    add_regression(axs[0][1], x, y)
    axs[0][1].set_title("CHM13 SVs")

    # Bottom left: alt contig SNV
    y = np.array([c.alt_snvs for c in counts])
    axs[1][0].scatter(x, y)
    add_regression(axs[1][0], x, y)
    axs[1][0].set_title("Alt contig SNVs")
    axs[1][0].set_ylabel("# variants")
    axs[1][0].set_xlabel("Distance to CHM13")

    # Bottom right: alt contig SV
    y = np.array([c.alt_svs for c in counts])
    axs[1][1].scatter(x, y)
    add_regression(axs[1][1], x, y)
    axs[1][1].set_title("Alt contig SVs")
    axs[1][1].set_xlabel("Distance to CHM13")

    fig.tight_layout()
    fig.savefig(f"{output_prefix}.n{n}.png")
    plt.close(fig)


if __name__ == "__main__":
    args = parse_args()

    dist_dict = load_distances(args.dist_matrix)

    with open(args.path_list) as f:
        paths = [line.strip() for line in f if line.strip()]

    # Store per-haplotype aggregated data
    hap_data = {n: {"dist": [], "counts": []} for n in range(1, 9)}

    for path in paths:
        logfile = os.path.join(args.log_dir, f"{path}.log")

        if not os.path.exists(logfile):
            raise FileNotFoundError(f"Missing logfile: {logfile}")

        if path not in dist_dict:
            raise KeyError(f"No CHM13 distance found for path: {path}")

        all_counts = parse_logfile(logfile)

        for n in range(1, 9):
            hap_data[n]["dist"].append(dist_dict[path])
            hap_data[n]["counts"].append(all_counts[n-1])

    # Generate plots
    for n in range(1, 9):
        plot_for_haplotype(
            n,
            hap_data[n]["dist"],
            hap_data[n]["counts"],
            args.output_prefix
        )