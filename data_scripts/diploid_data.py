#!/usr/bin/env python3
"""Collect data about cenhap alignment attempts.

== Output ==

For each diploid sample run, collect
- Chromosome
- Sample name
- Truth cenhaps
- Guessed cenhaps
- Stats (identity, correctness, runtime, memory usage) for each alignment run:
    - Minimap2 CHM13
    - Giraffe neighbor
    - Giraffe sampled
    - Giraffe native hap

These are output as a TSV to standard output.

== Inputs ==

General inputs, processed before specific haplotypes:
- Sample list (--sample-list).
    A CSV (with header) whose first column is a list of samples
- Truth cenhaps (<--cenhap-dir>/<chrom>/<chrom>.cenhap_predictions.tsv).
    Expects columns for path name and then cenhap, in that order.
    Also expects a header line but it will be skipped.

Sample-specific inputs:
- Haplotype sampling logs (<--log-dir>/<chrom>.<sample>.guess.real.log).
    Uses lines like "Selected haplotype <name> with score <score>"
- Alignment stats (<--aln-dir>/<chrom>.<path>.stats.log).
    A stats file with lines
        [prefix]: identity [0-1] / correctness [%] / runtime [sec] / memory [GB]

== How I ran this ==

From within /private/home/fokamoto/centromere-haplotype-sampling-pipeline

./data_scripts/diploid_data.py \
    -c input_data \
    -s input_data/aws_file_locations.csv \
    -l /private/groups/patenlab/fokamoto/centrolign/graph/diploid \
    -a /private/groups/patenlab/fokamoto/centrolign/alignments/diploid > diploid_data.tsv
"""

import argparse # Command-line argument parsing
from collections import Counter # Easier counting of cenhap matches
import os # Filesystem interactions
import sys # stderr
from typing import Dict, List, Tuple # Type hinting

# Parts of file names (see file docstring)
CENHAP_SUFFIX = 'cenhap_predictions.tsv'
GUESS_SUFFIX = 'guess.real.log'
STATS_SUFFIX = 'stats.log'

REFS = {'CHM13'      : 'CHM13', 
        'Neighbor'   : 'neighbor', 
        'Sampled'    : 'sampled', 
        'Native hap' : 'own_hap'}
"""Match up column titles to filenames"""
# All combinations for diploid alignments
ALN_COMBOS = [('CHM13', 'real', 'minimap2'),
              ('Neighbor', 'real', 'giraffe'),
              ('Sampled', 'real', 'giraffe'),
              ('Native hap', 'real', 'giraffe')]

# A full list of chromosomes to search for
CHROMOSOMES = [f'chr{i}' for i in range(23)] + ['chrX']

Diplotype = Tuple[str, str]
"""A sorted pair of cenhaps assigned to a sample."""

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Get diploid typing data.')
    parser.add_argument('-c', '--cenhap-dir', required=True,
                        help='Directory with TSVs of truth cenhap assignments')
    parser.add_argument('-s', '--sample-list', required=True,
                        help='CSVs with sample names')
    parser.add_argument('-l', '--log-dir', required=True,
                        help='Directory with haplotype sampling logs')
    parser.add_argument('-a', '--aln-dir', required=True,
                        help='Directory with read alignment stat files')
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

def read_sample_list(sample_list: str) -> List[str]:
    """Read samples from the first column of a CSV."""
    samples = []
    with open(sample_list) as file:
        file.readline()  # Skip header
        for line in file:
            samples.append(line.strip().split(',')[0])
    return samples

def read_all_cenhap_tables(cenhap_dir: str) -> Dict[str, Dict[str, str]]:
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

def construct_diplotype(cenhap_1: str, cenhap_2: str) -> Diplotype:
    """Normalize a diplotype by sorting the two cenhaps."""
    if cenhap_1 < cenhap_2:
        return cenhap_1, cenhap_2
    else:
        return cenhap_2, cenhap_1

def extract_guess(log_file: str) -> Tuple[str, str]:
    """Read haplotype pair from file."""
    with open(log_file) as file:
        for line in file:
            parts = line.strip().split()
            if line.startswith('Best guess'):
                return construct_diplotype(parts[9], parts[11])
    raise ValueError('No cenhap guess found in log file')

def count_matches(guess: Diplotype, truth: Diplotype) -> int:
    """Count how many cenhaps match between two diplotypes."""
    guess_counter = Counter(guess)
    truth_counter = Counter(truth)
    return sum(min(truth_counter[cenhap], count)
               for cenhap, count in guess_counter.items())

def get_aln_stats(log_file: str) -> Dict[str, List[float]]:
    """Read alignment statistics from stats log file.
    
    Uses lines like
        [prefix]: identity [0-1] / correctness [%] / runtime [sec] / memory [GB]
    
    Returns {prefix : [identity, correctness, runtime, memory]}
    """

    aln_stats = dict()
    with open(log_file) as file:
        for line in file:
            if not 'identity' in line:
                # Not alignment stats
                continue
            parts = line.strip().split()

            prefix = parts[0].rstrip(':')
            identity = float(parts[2])
            if parts[5] == 'None':
                correctness = None
            else:
                correctness = float(parts[5])
            runtime = float(parts[8])
            memory = float(parts[11])

            aln_stats[prefix] = [identity, correctness, runtime, memory]

    return aln_stats

if __name__ == '__main__':
    args = parse_args()
    cenhap_tables = read_all_cenhap_tables(args.cenhap_dir)
    samples = sorted(read_sample_list(args.sample_list))

    column_titles = ['Chromosome', 'Sample', 
                     'Truth cenhaps', 'Guessed cenhaps']
    for (ref, realness, tool) in ALN_COMBOS:
        aln_group = f'{ref} {tool} {realness}'
        column_titles += [f'{aln_group} identity', f'{aln_group} correctness',
                          f'{aln_group} runtime', f'{aln_group} memory']
    print('\t'.join(column_titles))

    for chrom in CHROMOSOMES:
        exact_matches = 0
        near_matches = 0
        total = 0
        
        for name in samples:
            items_to_write = [chrom, name]

            try:
                true_cenhaps = construct_diplotype(
                    cenhap_tables[chrom][f'{name}.1'],
                    cenhap_tables[chrom][f'{name}.2'])
                
                guess_log = os.path.join(args.log_dir, 
                                        f'{chrom}.{name}.{GUESS_SUFFIX}')
                guessed_cenhaps = extract_guess(guess_log)
                items_to_write += [','.join(true_cenhaps),
                                    ','.join(guessed_cenhaps)]
                
                match_count = count_matches(guessed_cenhaps, true_cenhaps)
                if match_count == 2:
                    exact_matches += 1
                elif match_count == 1:
                    near_matches += 1
                total += 1
            except:
                # Failed to get cenhap data
                items_to_write += ['None', 'None']
            
            # Dump in alignment stats as well
            stats_file = os.path.join(
                args.aln_dir, f'{chrom}.{name}.{STATS_SUFFIX}')
            if os.path.exists(stats_file):
                aln_stats = get_aln_stats(os.path.join(
                    args.aln_dir, f'{chrom}.{name}.{STATS_SUFFIX}'))
                if len(ALN_COMBOS) != len(aln_stats):
                    continue
                for (ref, realness, tool) in ALN_COMBOS:
                    prefix = f'{chrom}.{name}.{REFS[ref]}.{realness}.{tool}'
                    items_to_write += aln_stats[prefix]
                
                print('\t'.join(str(item) for item in items_to_write))

        if chrom in cenhap_tables:
            exact_acc = 100 * exact_matches / total
            near_acc = 100 * (exact_matches + near_matches) / total
            print(f'{chrom} accuracy: {near_acc:.2f}% half, {exact_acc:.2f}% exact',
                file=sys.stderr)