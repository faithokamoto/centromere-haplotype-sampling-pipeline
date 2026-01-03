#!/usr/bin/env python3
"""Plot identity and accuracy across n haplotypes sampled.

Reads its data from a logfile.
"""

import argparse
from dataclasses import dataclass
from typing import List

import matplotlib.pyplot as plt

@dataclass
class AccuracyStats:
    """Accuracy statistics."""
    precision: float
    recall: float
    recall_excl_svs: float

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot identity and accuracy from logfile."
    )
    parser.add_argument(
        "--name", type=str, help="Sample/haplotype name."
    )
    parser.add_argument(
        "--logfile",
        type=str,
        help="Path to the logfile containing identity and accuracy stats."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        help="Path to save the output plot."
    )
    return parser.parse_args()

def get_identity_stats(logfile: str) -> List[int]:
    """Get identity stats from logfile."""
    avg_identities = []
    with open(logfile, 'r') as f:
        for line in f:
            if "avg identity=" in line:
                parts = line.strip().split('avg identity=')
                identity = float(parts[1].split()[0])
                avg_identities.append(identity)
    return avg_identities

def get_snv_stats(logfile: str) -> List[AccuracyStats]:
    """Get SNV stats from logfile."""
    snv_stats = []
    with open(logfile, 'r') as f:
        for line in f:
            if line.startswith("Precision:"):
                precision = float(line.strip().split()[1])
            elif line.startswith("Recall:"):
                recall = float(line.strip().split()[1])
            elif line.startswith("Recall excluding FNs in SVs:"):
                recall_excl_svs = float(line.strip().split()[5])
                snv_stats.append(AccuracyStats(precision, recall, recall_excl_svs))
    return snv_stats

def plot_identity_and_accuracy(id_stats: List[int], 
                               snv_stats: List[AccuracyStats], 
                               name: str,
                               output_file: str) -> None:
    """Plot identity and accuracy from logfile."""
    title = f'{name} stats across number of haplotypes sampled'

    ns = list(range(1, len(id_stats) + 1))

    fig, ax1 = plt.subplots(figsize=(8, 5))

    color = 'tab:blue'
    ax1.set_xlabel('Number of Haplotype Samples (n)')
    ax1.set_ylabel('Average Identity', color=color)
    ax1.plot(ns, id_stats, 'o-', color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    if len(snv_stats) == 0:
        print("No SNV stats found in logfile; skipping SNV accuracy plot.")
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(output_file)
        return

    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('SNV Accuracy', color=color)
    precisions = [stat.precision for stat in snv_stats]
    recalls = [stat.recall for stat in snv_stats]
    recalls_excl_svs = [stat.recall_excl_svs for stat in snv_stats]
    ax2.plot(ns, precisions, 's--', label='Precision', color='tab:orange')
    ax2.plot(ns, recalls, 'd--', label='Recall', color='tab:red')
    ax2.plot(ns, recalls_excl_svs, 'x--', label='Recall Excl. SVs', color='tab:pink')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.legend(loc='upper left')

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_file)

if __name__ == "__main__":
    args = parse_args()
    identity_stats = get_identity_stats(args.logfile)
    snv_stats = get_snv_stats(args.logfile)
    plot_identity_and_accuracy(identity_stats, snv_stats, 
                               args.name, args.output_file)