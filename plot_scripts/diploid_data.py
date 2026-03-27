#!/usr/bin/env python3
"""Collect data about cenhap alignment attempts.

== Output ==

For each diploid sample run, collect
- Chromosome
- Sample name
- Truth cenhaps
- Guessed cenhaps

These are output as a TSV to standard output.

== Inputs ==

General inputs, processed before specific haplotypes:
- Truth cenhaps (<--cenhap-dir>/<chrom>/<chrom>.cenhap_predictions.tsv).
    Expects columns for path name and then cenhap, in that order.
    Also expects a header line but it will be skipped.

Sample-specific inputs:
- Haplotype sampling logs (<--log-dir>/<chrom>.<sample>.guess.real.log).
    Expects a line output by `./guess_n_and_cenhap.py`.
    "Best guess: use <n> haplotypes & sample cenhap = <cenhap>".
"""

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Dict # Type hinting

# Parts of file names (see file docstring)
CENHAP_SUFFIX = 'cenhap_predictions.tsv'
GUESS_SUFFIX = 'guess.real.log'

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Get diploid typing data.')
    parser.add_argument('-c', '--cenhap-dir', required=True,
                        help='Directory with TSVs of truth cenhap assignments')
    parser.add_argument('-l', '--log-dir',
                        help='Directory with haplotype sampling logs')
    return parser.parse_args()

def read_chrom_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read a cenhap assignment table.
    
    Assuming the first line is a header,
    reads haplotype,cenhap lines into 
    {sample : "cenhap1,cenhap2"} dictionary.
    """
    cenhap_table = dict()
    all_samples = set()
    # Collect all hap -> cenhap assignments
    with open(cenhap_file) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            hap_name = parts[0]
            sample_name = hap_name.split('.')[0]
            cenhap_table[hap_name] = parts[1]
            all_samples.add(sample_name)

    # Get diploid cenhaps when possible
    sample_to_cenhaps = dict()
    for sample in all_samples:
        if f'{sample}.1' in cenhap_table and f'{sample}.2' in cenhap_table:
            cenhap1 = cenhap_table[f'{sample}.1']
            cenhap2 = cenhap_table[f'{sample}.2']
            sample_to_cenhaps[sample] = f'{cenhap1},{cenhap2}'
    return sample_to_cenhaps

def read_all_cenhap_tables(cenhap_dir: str) -> Dict[str, Dict[str, str]]:
    """Reads all cenhap tables in a directory.
    
    Finds all files matching
        <--cenhap-dir>/<chrom>/<chrom>.cenhap_predictions.tsv
    and reads them into a {chrom : table} dictionary
    """
    cenhap_tables = dict()
    for item in os.listdir(cenhap_dir):
        if 'chr' in item:
            cenhap_tables[item] = read_chrom_cenhap_table(
                os.path.join(cenhap_dir, item, f'{item}.{CENHAP_SUFFIX}'))
    return cenhap_tables

def extract_guess(log_path):
    """Extract cenhaps from 'Best guess <etc> 1 / 2'"""
    with open(log_path) as f:
        for line in f:
            if line.startswith('Best guess'):
                return f'{line.split()[-3]},{line.split()[-1]}'

if __name__ == '__main__':
    args = parse_args()
    sample_to_true = read_all_cenhap_tables(args.cenhap_dir)

    print('\t'.join(['Chromosome', 'Sample', 
                     'Truth cenhaps', 'Guessed cenhaps']))
    for chrom in sample_to_true.keys():
        for sample, true_cenhaps in sample_to_true[chrom].items():
            guess_log = os.path.join(args.log_dir, 
                                     f'{chrom}.{sample}.{GUESS_SUFFIX}')
            if os.path.exists(guess_log):
                guessed_cenhaps = extract_guess(guess_log)
                print('\t'.join([chrom, sample, true_cenhaps, guessed_cenhaps]))