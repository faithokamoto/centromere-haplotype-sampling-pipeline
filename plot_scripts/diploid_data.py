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
"""

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Dict, Set # Type hinting

# Parts of file names (see file docstring)
CENHAP_SUFFIX = 'cenhap_predictions.tsv'
GUESS_SUFFIX = 'guess.real.log'
DROP_THRESHOLD = 2000

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
    {haplotype : cenhap} dictionary.
    """
    cenhap_table = dict()
    # Collect all hap -> cenhap assignments
    with open(cenhap_file) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            hap_name = parts[0]
            cenhap_table[hap_name] = parts[1]
    return cenhap_table

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

def extract_guess(cenhap_tables: Dict[str, Dict[str, str]], chrom: str,
                  log_file: str) -> str:
    """Guess haplotype pair."""
    guessed_cenhaps = set()
    top_score = None
    with open(log_file) as file:
        for line in file:
            parts = line.strip().split()
            if line.startswith('Selected haplotype'):
                hap_name = parts[2]
                score = float(parts[5])
                if top_score is None:
                    top_score = score
                
                if top_score - score > DROP_THRESHOLD:
                    break

                if len(guessed_cenhaps) < 2:
                    if hap_name in cenhap_tables[chrom]:
                        this_cenhap = cenhap_tables[chrom][hap_name]
                        if not this_cenhap in guessed_cenhaps:
                            guessed_cenhaps.add(this_cenhap)
    
    guessed_cenhaps = sorted(guessed_cenhaps)
    if len(guessed_cenhaps) == 1:
        return f'{guessed_cenhaps[0]},{guessed_cenhaps[0]}'
    elif len(guessed_cenhaps) == 2:
        return f'{guessed_cenhaps[0]},{guessed_cenhaps[1]}'
    else:
        return None

if __name__ == '__main__':
    args = parse_args()
    cenhap_tables = read_all_cenhap_tables(args.cenhap_dir)

    print('\t'.join(['Chromosome', 'Sample', 
                     'Truth cenhaps', 'Guessed cenhaps']))
    for chrom in cenhap_tables.keys():
        cur_table = cenhap_tables[chrom]
        all_samples = {hap.split('.')[0] for hap in cur_table.keys()}
        sample_to_true = dict()
        for sample in all_samples:
            if f'{sample}.1' in cur_table and f'{sample}.2' in cur_table:
                cenhap1 = cur_table[f'{sample}.1']
                cenhap2 = cur_table[f'{sample}.2']
                truth = sorted([cenhap1, cenhap2])
                sample_to_true[sample] = f'{truth[0]},{truth[1]}'
        
        for sample, true_cenhaps in sample_to_true.items():
            guess_log = os.path.join(args.log_dir, 
                                     f'{chrom}.{sample}.{GUESS_SUFFIX}')
            if os.path.exists(guess_log):
                guessed_cenhaps = extract_guess(cenhap_tables, chrom, guess_log)
                print('\t'.join([chrom, sample, true_cenhaps, guessed_cenhaps]))