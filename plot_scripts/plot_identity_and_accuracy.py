#!/usr/bin/env python3
"""Plot stats for CHM13 & alt contigs across n haplotypes sampled.

Assumes the logfile is for a single sample and has the format:

...
Sampling 1 haplotypes
...
1791 SNVs and 446 SVs on CHM13
518 SNVs and 214 SVs on alt contigs
...
Overall:
True Positives (TP): 188
False Positives (FP): 139
Precision: 0.5749

On alt contigs:
True Positives (TP): 188
False Positives (FP): 139
Precision: 0.5749
...
Sampling 2 haplotypes
...

n=1: avg identity=0.991656
n=2: avg identity=0.994935
  Significant jump detected at n=2!
n=3: avg identity=0.995063
n=4: avg identity=0.994998
n=5: avg identity=0.995170
n=6: avg identity=0.995706
n=7: avg identity=0.995810
n=8: avg identity=0.995826
Sample 2 haplotypes for HG01891_pat
"""

import argparse # Command line interface
from dataclasses import dataclass # Simple struct
from typing import List, Tuple # Type hinting

import matplotlib.pyplot as plt # Basic plotting

@dataclass
class RunStats:
    """Variant statistics for a haplotype sampling run"""
    # Count how many variants were called
    chm13_snvs: int
    chm13_svs : int
    alt_snvs: int
    alt_svs: int

    # Count TP/FP for ones with a truth set
    chm13_tp: int = 0
    chm13_fp: int = 0
    alt_tp: int = 0
    alt_fp: int = 0

    # Overall read identity
    avg_identity: float = 0

    def set_tp_fp(self, chm13_tp: int, chm13_fp: int, 
                  alt_tp: int, alt_fp: int) -> None:
        """Set all the TP/FP fields"""
        self.chm13_tp = chm13_tp
        self.chm13_fp = chm13_fp
        self.alt_tp = alt_tp
        self.alt_fp = alt_fp

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
        "--output-file",
        type=str,
        help="Path to save the output plot."
    )
    return parser.parse_args()

def get_snv_stats(logfile: str) -> Tuple[List[RunStats], int]:
    """Get stats from logfile.
    
    Creates a RunStats for each haplotype sampling run,
    as well as grabbing the optimal n from the end of the file.
    Assumes a very strict format; see file docstring.
    """

    stats = []
    chosen_n = None
    with open(logfile, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].strip().endswith('SVs on CHM13'):
                # Detect when we're reporing the number of variants called
                chm13_snvs = int(lines[i].split()[0])
                chm13_svs = int(lines[i].split()[3])
                alt_snvs = int(lines[i+1].split()[0])
                alt_svs = int(lines[i+1].split()[3])
                stats.append(RunStats(chm13_snvs, chm13_svs, alt_snvs, alt_svs))
            elif lines[i].startswith("Overall"):
                # Overal stats header means individual stats lines are nearby
                overall_tp = int(lines[i+1].split()[3])
                overall_fp = int(lines[i+2].split()[3])
                alt_tp = int(lines[i+6].split()[3])
                alt_fp = int(lines[i+7].split()[3])
                chm13_tp = overall_tp - alt_tp
                chm13_fp = overall_fp - alt_fp
                stats[-1].set_tp_fp(chm13_tp, chm13_fp, alt_tp, alt_fp)
            elif 'avg identity=' in lines[i]:
                parts = lines[i].strip().split()
                num_hap = int(parts[0].split('=')[1].strip(':'))
                identity = float(parts[2].split('=')[1])
                stats[num_hap-1].avg_identity = identity
            elif 'haplotypes for' in lines[i]:
                parts = lines[i].split()
                chosen_n = int(parts[1])
    return stats, chosen_n

def plot_identity_and_accuracy(stats: List[RunStats], chosen_n: int,
                               name: str,
                               output_file: str) -> None:
    """Plot identity and accuracy from logfile."""
    ns = list(range(1, len(stats) + 1))

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

    # Top left: CHM13 SNV TP/FP
    axs[0][0].plot(ns, [s.chm13_fp for s in stats], 'o-', 
                   color='#E69F00', label='FP')
    axs[0][0].plot(ns, [s.chm13_tp for s in stats], 'o-', 
                   color='#56B4E9', label='TP')
    axs[0][0].set_ylabel('# SNVs')
    axs[0][0].set_title(f'{name} vs. CHM13')
    axs[0][0].legend(loc='upper right')

    # Top right: alt contig SNV TP/FP
    axs[0][1].plot(ns, [s.alt_fp for s in stats], 'o-', 
                   color='#E69F00', label='FP')
    axs[0][1].plot(ns, [s.alt_tp for s in stats], 'o-', 
                   color='#56B4E9', label='TP')
    axs[0][1].set_title(f'{name} vs. alts')
    axs[0][1].legend(loc='upper right')

    # Bottom left: raw variant counts
    axs[1][0].plot(ns, [s.chm13_snvs for s in stats], 'o-', 
                   color='#009E73', label='CHM13 SNVs')
    axs[1][0].plot(ns, [s.chm13_svs for s in stats], 'o-', 
                   color='#0072B2', label='CHM13 SVs')
    axs[1][0].plot(ns, [s.alt_snvs for s in stats], 'o-', 
                   color='#D55E00', label='alt contig SNVs')
    axs[1][0].plot(ns, [s.alt_svs for s in stats], 'o-', 
                   color='#CC79A7', label='alt contig SVs')
    axs[1][0].set_xlabel('# haplotypes sampled')
    axs[1][0].set_ylabel('# variants')
    axs[1][0].set_title('Variant count')
    axs[1][0].legend(loc='upper right')

    # Bottom right: average read identity
    axs[1][1].plot(ns, [s.avg_identity for s in stats], 'o-',
                   color='black')
    axs[1][1].plot(chosen_n, stats[chosen_n-1].avg_identity, '*',
                   color='#F0E442')
    axs[1][1].set_xlabel('# haplotypes sampled')
    axs[1][1].set_title('Average read identity')

    fig.tight_layout()
    fig.savefig(output_file)

if __name__ == "__main__":
    args = parse_args()
    snv_stats, chosen_n = get_snv_stats(args.logfile)
    plot_identity_and_accuracy(snv_stats, chosen_n, args.name, args.output_file)