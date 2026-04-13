#!/usr/bin/env python3
"""Make Panel C of --absent-score supplementary figure.

Go through haplotype sampling log files in a directory,
and plot # of times HG00738.1 tricked the typer.

Uses lines like
    Selected haplotype <haplotype> with score <score>
from files <--log-dir>/chr4.<haplotype>.<absent-score>.log

./param_tests/incorrect_vs_param.py \
    -c /private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/chr4/chr4.cenhap_predictions.tsv \
    -l log/typing_tests -o plot_outputs/incorrect_vs_param.png -n HG00738.1
"""

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Dict # Type hinting

import matplotlib.pyplot as plt # Basic plotting

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cenhap-table', required=True,
                        help='TSV of truth cenhap assignments')
    parser.add_argument('-l', '--log-dir', required=True,
                        help='Directory with log files')
    parser.add_argument('-o', '--output-file', required=True,
                        help='Graph file to write to')
    parser.add_argument('-n', '--haplotype-name', required=True,
                        help='Haplotype to call out')
    return parser.parse_args()

def read_chrom_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read a cenhap assignment table."""
    cenhap_table = dict()
    with open(cenhap_file) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            cenhap_table[parts[0]] = parts[1]
    return cenhap_table

def parse_filename(filename: str):
    """
    Parse filename of form:
    chr4.<haplotype>.<##>.log

    Returns:
        haplotype (str), absent_score (float)
    """
    parts = filename.split('.')
    # parts = ["chr4", sample_id, hap_num, "##", ".log"]
    score_part = parts[3]
    haplotype = '.'.join(parts[1:3])

    return haplotype, (int(score_part) / 100.0)

def parse_log_file(filepath: str):
    """
    Extract selected haplotype from log file.
        Selected haplotype <other hap> with score <score>
    """
    with open(filepath) as file:
        for line in file:
            if 'Selected' in line:
                return line.strip().split()[2]

if __name__ == '__main__':
    args = parse_args()
    cenhap_table = read_chrom_cenhap_table(args.cenhap_table)

    results = dict()

    for filename in os.listdir(args.log_dir):
        if not filename.endswith('log'):
            continue

        typing_hap, score = parse_filename(filename)
        if not score in results:
            results[score] = 0

        selected_hap = parse_log_file(os.path.join(args.log_dir, filename))
        true_cenhap = cenhap_table[typing_hap]
        selected_cenhap = cenhap_table[selected_hap]

        typing_success = (true_cenhap == selected_cenhap)

        if selected_hap == args.haplotype_name and not typing_success:
            results[score] += 1

    scores = sorted(results.keys())
    special_errors = [results[s] for s in scores]

    plt.figure()
    plt.scatter(scores, special_errors)
    plt.xlabel('--absent-score')
    plt.ylabel(f'# wrong {args.haplotype_name}')

    plt.savefig(args.output_file, bbox_inches='tight')
