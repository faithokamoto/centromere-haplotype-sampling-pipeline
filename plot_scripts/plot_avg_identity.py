#!/usr/bin/env python3
"""Plot average identity in a few conditions.

Compare CHM13, optimal haplotype sampled graph, and native haplotype.
"""

import csv
from typing import Dict

import matplotlib.pyplot as plt

# Where the alignment TSVs live
PROJ_DIR = '/private/groups/patenlab/fokamoto/centrolign'
LINEAR_DIR = f'{PROJ_DIR}/alignments/linear_refs/'
SAMPLED_DIR = f'{PROJ_DIR}/alignments/leave_one_out/'

# I/O file locations
INPUT_FILE = 'input_data/test_samples.txt'
OUTPUT_FILE = 'plot_outputs/avg_identities.png'

LABELS = ['Native Haplotype', 'Optimal # Sampled', 
          'Matrix Neighbor', '1 sampled', 'CHM13']
"""Conditions to plot, from bottom to top"""

def get_identity_tsv_avg(tsv_file: str) -> float:
    """Average read identity from a TSV file."""
    identities = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            identities.append(float(row['identity']))
    return sum(identities) / len(identities)

def get_sample_avgs(sample_name: str, optimal_hap: int) -> Dict[str, float]:
    """Get average identities for a given sample."""
    avgs = dict()
    filenames = [f'{LINEAR_DIR}/{sample_name}.own_hap.real.giraffe.tsv',
                 f'{SAMPLED_DIR}/real_{sample_name}.{optimal_hap}haps.tsv',
                 f'{LINEAR_DIR}/{sample_name}.neighbor.real.giraffe.tsv',
                 f'{SAMPLED_DIR}/real_{sample_name}.1haps.tsv',
                 f'{LINEAR_DIR}/{sample_name}.chm13.real.giraffe.tsv']
    avgs[LABELS[0]] = get_identity_tsv_avg(filenames[0])

    if optimal_hap > 0:
        # Only use haplotype sampling for non-hopeless samples
        avgs[LABELS[1]] = get_identity_tsv_avg(filenames[1])
    else:
        # Otherwise, use 1 sampled haplotype
        avgs[LABELS[1]] = get_identity_tsv_avg(filenames[3])
    
    avgs[LABELS[2]] = get_identity_tsv_avg(filenames[2])
    avgs[LABELS[3]] = get_identity_tsv_avg(filenames[3])
    avgs[LABELS[4]] = get_identity_tsv_avg(filenames[4])
    return avgs

def get_all_sample_avgs(sample_file: str) -> Dict[str, Dict[str, float]]:
    """Get average identities for all samples.

    Parameters
    ----------
    sample_file : str
        Path to the sample file.

    Returns
    -------
    all_avgs : Dict[str, Dict[str, float]]
        {{ sample name : { condition : avg identity } }
    """
    all_avgs = dict()
    with open(sample_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            all_avgs[parts[2]] = get_sample_avgs(parts[2], int(parts[3]))
    return all_avgs

def plot_avgs(all_avgs: Dict[str, Dict[str, float]], output_file: str) -> None:
    """Plot the average identities for different conditions.

    Parameters
    ----------
    all_avgs : Dict[str, Dict[str, float]]
        {{ sample name : { condition : avg identity } }
    output_file : str
        Path to save the output plot.
    """

    plt.figure(figsize=(4, 3))
    
    for avg in all_avgs.values():
        x = [avg[condition] for condition in LABELS if condition in avg]
        y = [idx for idx, condition in enumerate(LABELS) if condition in avg]
        plt.plot(x, y, 'o-', color='gray', linewidth=0.5)
        plt.plot(x, y, 'o', color='black')

    plt.yticks(range(len(LABELS)), LABELS)
    plt.xlabel('Average Identity')
    plt.tight_layout()
    plt.savefig(output_file)

if __name__ == '__main__':
    all_avgs = get_all_sample_avgs(INPUT_FILE)
    plot_avgs(all_avgs, OUTPUT_FILE)