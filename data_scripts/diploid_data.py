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
    Uses lines like "Selected haplotype <name> with score <score>"

== How I ran this ==

From within /private/home/fokamoto/centromere-haplotype-sampling-pipeline

./data_scripts/diploid_data.py \
    -c input_data \
    -l /private/groups/patenlab/fokamoto/centrolign/graph/diploid > diploid_data.tsv
"""

import argparse # Command-line argument parsing
from collections import Counter # Easier counting of cenhap matches
import os # Filesystem interactions
import sys # stderr
from typing import Dict, Tuple # Type hinting

# Parts of file names (see file docstring)
CENHAP_SUFFIX = 'cenhap_predictions.tsv'
GUESS_SUFFIX = 'guess.real.log'

Diplotype = Tuple[str, str]
"""A sorted pair of cenhaps assigned to a sample."""

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Get diploid typing data.')
    parser.add_argument('-c', '--cenhap-dir', required=True,
                        help='Directory with TSVs of truth cenhap assignments')
    parser.add_argument('-l', '--log-dir', required=True,
                        help='Directory with haplotype sampling logs')
    return parser.parse_args()

def read_chrom_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read a cenhap assignment table.
    
    Assuming the first line is a header,
    reads haplotype,cenhap lines into 
    {haplotype : cenhap} dictionary.
    """
    cenhap_table = dict()
    with open(cenhap_file) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            cenhap_table[parts[0]] = parts[1]
    return cenhap_table

def read_all_cenhap_tables(cenhap_dir: str) -> Dict[str, Diplotype]:
    """Reads all cenhap tables in a directory.
    
    Finds all files matching
        <--cenhap-dir>/<chrom>.cenhap_predictions.tsv
    and reads them into a {chrom : table} dictionary
    """
    cenhap_tables = dict()
    for item in os.listdir(cenhap_dir):
        if CENHAP_SUFFIX in item:
            chrom = item.split('.')[0]
            cenhap_tables[chrom] = read_chrom_cenhap_table(
                os.path.join(cenhap_dir, item))
    return cenhap_tables

def extract_guess(log_file: str) -> Tuple[str]:
    """Read haplotype pair from file."""
    with open(log_file) as file:
        for line in file:
            parts = line.strip().split()
            if line.startswith('Best guess'):
                guess_1 = parts[9]
                guess_2 = parts[11]
                return tuple(sorted([guess_1, guess_2]))
            
def find_all_diplotypes(cur_table: Dict[str, str]) -> Dict[str, Diplotype]:
    """Convert a {haplotype : cenhap} table to a {sample : cenhaps} table."""
    all_samples = {hap.split('.')[0] for hap in cur_table.keys()}
    sample_to_true = dict()
    for sample in all_samples:
        if f'{sample}.1' in cur_table and f'{sample}.2' in cur_table:
            cenhap1 = cur_table[f'{sample}.1']
            cenhap2 = cur_table[f'{sample}.2']
            sample_to_true[sample] = tuple(sorted([cenhap1, cenhap2]))
    return sample_to_true

def count_matches(guess: Diplotype, truth: Diplotype) -> int:
    """Count how many cenhaps match between two diplotypes."""
    guess_counter = Counter(guess)
    truth_counter = Counter(truth)
    return sum(min(truth_counter[cenhap], count)
               for cenhap, count in guess_counter.items())

if __name__ == '__main__':
    args = parse_args()
    cenhap_tables = read_all_cenhap_tables(args.cenhap_dir)

    print('\t'.join(['Chromosome', 'Sample', 
                     'Truth cenhaps', 'Guessed cenhaps']))
    for chrom in cenhap_tables.keys():
        sample_to_true = find_all_diplotypes(cenhap_tables[chrom])
        exact_matches = 0
        near_matches = 0
        total = 0

        for sample, true_cenhaps in sample_to_true.items():
            guess_log = os.path.join(args.log_dir, 
                                     f'{chrom}.{sample}.{GUESS_SUFFIX}')
            if os.path.exists(guess_log):
                guessed_cenhaps = extract_guess(guess_log)
                print('\t'.join([chrom, sample, 
                                ','.join(true_cenhaps),
                                ','.join(guessed_cenhaps)]))
                
                match_count = count_matches(guessed_cenhaps, true_cenhaps)
                if match_count == 2:
                    exact_matches += 1
                elif match_count == 1:
                    near_matches += 1
                total += 1
        
        exact_acc = 100 * exact_matches / total
        near_acc = 100 * (exact_matches + near_matches) / total
        print(f'{chrom} accuracy: {near_acc:.2f}% half, {exact_acc:.2f}% exact',
              file=sys.stderr)