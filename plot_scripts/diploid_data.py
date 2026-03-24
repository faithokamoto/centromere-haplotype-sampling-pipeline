#!/usr/bin/env python3

import argparse # Command-line argument parsing
import os # Filesystem interactions

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Get diploid typing data.')
    parser.add_argument('-c', '--cenhap-table',
                        help='TSV: haplotype (<sample>.1/.2) and cenhap')
    parser.add_argument('-l', '--log-dir',
                        help='Directory with <sample>.chr12.guess.real.log')
    return parser.parse_args()

def load_cenhaps(path: str) -> dict:
    """Look up truth cenhaps from table, e.g. 1,2"""
    sample_to_cenhaps = {}

    with open(path) as f:
        f.readline() # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            hap, cenhap = line.split()
            sample, _= hap.rsplit('.', 1)
            sample_to_cenhaps.setdefault(sample, []).append(cenhap)

    # Convert formats
    return {sample: ','.join(cenhaps) 
            for sample, cenhaps in sample_to_cenhaps.items()}

def extract_guess(log_path):
    """Extract cenhaps from 'Best guess <etc> 1 / 2'"""
    with open(log_path) as f:
        for line in f:
            if line.startswith('Best guess'):
                return f'{line.split()[-3]},{line.split()[-1]}'

if __name__ == '__main__':
    args = parse_args()
    sample_to_true = load_cenhaps(args.cenhap_table)

    print('\t'.join(['sample', 'guessed_cenhaps', 'true_cenhaps']))
    for sample, true_cenhaps in sorted(sample_to_true.items()):
        log_file = os.path.join(args.log_dir, f'{sample}.chr12.guess.real.log')
        if not os.path.exists(log_file):
            continue
        guessed = extract_guess(log_file)
        print("\t".join([sample, guessed, true_cenhaps]))