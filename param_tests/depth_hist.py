#!/usr/bin/env python3

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Tuple # Type hinting

import matplotlib.pyplot as plt # Basic plotting

GUESS_SUFFIX = '.guess.log'

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log-dir', required=True,
                        help='Directory with log files')
    parser.add_argument('-o', '--output-file', required=True,
                        help='Graph file to write to')
    parser.add_argument('-n', '--haplotype-name', required=True,
                        help='Haplotype to call out')
    return parser.parse_args()

def read_top_hap_stats(logfile: str) -> Tuple[str, float]:
    """Read the top haplotype & its private depth from logfile.
    
    Looks for a line like
        <hap name> has avg depth <depth>
    and returns (name, depth)
    """
    with open(logfile) as file:
        for line in file:
            if 'has avg depth' in line:
                parts = line.strip().split()
                top_name = parts[0]
                top_depth = float(parts[4])
                return top_name, top_depth
    return None, None

if __name__ == '__main__':
    args = parse_args()

    special_hap_depths = []
    other_hap_depths = []

    for item in os.listdir(args.log_dir):
        if item.endswith(GUESS_SUFFIX):
            hap, depth = read_top_hap_stats(os.path.join(args.log_dir, item))
            if hap is None:
                continue
            if hap == args.haplotype_name:
                special_hap_depths.append(depth)
            else:
                other_hap_depths.append(depth)

    # Plot grouped histogram
    plt.figure(figsize=(8, 6))
    plt.hist(
        [special_hap_depths, other_hap_depths],
        bins=[val / 2 for val in range(0, 61)],  # 0–30 inclusive
        label=[args.haplotype_name, 'Other'],
        alpha=0.7,
        edgecolor='black'
    )

    plt.xlabel('Depth on top haplotype')
    plt.ylabel('Haplotype count')
    plt.legend()

    plt.tight_layout()
    plt.savefig(args.output_file, dpi=150)
    plt.close()
