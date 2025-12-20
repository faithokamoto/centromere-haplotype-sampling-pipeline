#!/usr/bin/env python3
"""Guess optimal *n* for how many haplotypes to sample.

Input a haplotype name and this script will attempt to guess
what cenhap it should be, based on the haplotyeps sampled at
the optimal *n*. If it is hopeless, it will say so.
"""

import argparse

from typing import Dict, List, Tuple

SAMPLE_TABLE = 'input_data/test_samples.txt'
"""Four-column CSV of (path name, version, haplotype name, optimal N)."""
CENHAP_TABLE = 'input_data/chr12_r2_QC_v2_centrolign_cenhap_assignments.csv'
LOG_FILE = 'log/test_output.txt'

SampledHaplo = List[Tuple[str, float]]
"""List of (haplotype, score) tuples."""

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Guess optimal number of haplotypes to sample")
    parser.add_argument("haplotype", type=str,
                        help="Name of the haplotype (e.g., 'HG02622_mat')")
    return parser.parse_args()

def read_metadata(sample_table_file: str, haplo_name: str) -> Tuple[str, int]:
    """Read the optimal *n* & path name for a given haplotype"""
    with open(sample_table_file, 'r') as file:
        for line in file:
            path_name, _, name, optimal_n = line.strip().split(',')
            if name == haplo_name:
                return path_name, int(optimal_n)
    raise ValueError(f"Haplotype {haplo_name} not found in sample table.")

def read_sampled_haplotypes(log_file: str, haplo_name: str, n_hap: int) -> SampledHaplo:
    """Read haplotypes sampled for a given input in a log file."""
    sampled_haplotypes = []
    # How far have we read so far?
    found_haplo_name = False
    found_sample_run = False
    with open(log_file, 'r') as file:
        for line in file:
            # Finished that run
            if found_sample_run and ('autoindex' in line or 'Sampling' in line):
                break
            # Start of this sample's log
            if 'Processing sample' in line and haplo_name in line:
                found_haplo_name = True
            # Start of this sample's target sampling run
            if found_haplo_name and f'Sampling {n_hap} haplotypes' in line:
                found_sample_run = True
            # A haplotype sampled in the target sampling run
            if found_sample_run and 'Selected haplotype' in line:
                # Selected haplotype <name> with score <score>
                haplo = line.split()[2]
                score = float(line.split()[5])
                sampled_haplotypes.append((haplo, score))
    return sampled_haplotypes

def read_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read the cenhap assignment table."""
    cenhap_table = {}
    with open(cenhap_file, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            haplo, cenhap = parts[0], parts[1]
            cenhap_table[haplo] = cenhap
    return cenhap_table

def guess_cenhap(sampled: SampledHaplo, cenhap_table: Dict[str, str]) -> str:
    """Guess the cenhap for the sampled haplotypes."""
    # Normalize scores so minimum is 1
    min_score = min(score for _, score in sampled)
    normalized = [(haplo, score - min_score + 1) for haplo, score in sampled]
    # Calculate weighted score for each cenhap type
    cenhap_scores = {}
    for haplo, score in normalized:
        cenhap = cenhap_table.get(haplo, None)
        print("\tHaplotype:", haplo, "Cenhap:", cenhap, "Score:", score + min_score - 1)
        if cenhap is not None:
            cenhap_scores[cenhap] = cenhap_scores.get(cenhap, 0) + score
    # Guess the cenhap with the highest score
    return max(cenhap_scores.items(), key=lambda x: x[1])[0]

if __name__ == '__main__':
    args = parse_args()
    haplo_name = args.haplotype

    path_name, optimal_n = read_metadata(SAMPLE_TABLE, haplo_name)
    if optimal_n == 0:
        print(f"Haplotype {haplo_name} has optimal n=0; cannot guess cenhap.")
        exit(0)
    sampled_haplotypes = read_sampled_haplotypes(LOG_FILE, haplo_name, optimal_n)
    cenhap_table = read_cenhap_table(CENHAP_TABLE)

    guessed_cenhap = guess_cenhap(sampled_haplotypes, cenhap_table)
    correct_cenhap = cenhap_table.get(path_name, 'Unknown')
    print(f"Guessed cenhap for {haplo_name} (optimal n={optimal_n}): "\
          f"{guessed_cenhap} (correct: {correct_cenhap})")