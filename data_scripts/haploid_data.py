#!/usr/bin/env python3
"""Collect data about cenhap alignment attempts.

== Output ==

For each haplotype run, collect
- Chromosome
- Haplotype name
- Truth cenhap
- Guessed cenhap
- Names of sampled haps (comma-separated)
- Scores of sampled haps (comma-separated)
- Distance to closest hap in graph
- Distance to closest sampled hap (within the ones we chose to align to)
- Stats (identity, correctness, runtime, memory usage) for each alignment run:
    - Minimap2 CHM13
    - Giraffe CHM13
    - Minimap2 neighbor
    - Giraffe neighbor
    - Giraffe sampled
    - Minimap2 native hap
    - Giraffe native hap

These are output as a TSV to standard output.

== Inputs ==

General inputs, processed before specific haplotypes:
- Truth cenhaps (<--cenhap-dir>/<chrom>.cenhap_predictions.tsv).
    Expects columns for path name and then cenhap, in that order.
    Also expects a header line but it will be skipped.
- Pair dists (<--dist-dir>/<chrom>_r2_QC_v2_centrolign_pairwise_distance.csv).
    Expects columns for hap1, hap2, and then float distance (0-1).
    Does not expect a header line, or even symmetric lines.

Sample-specific inputs:
- Haplotype sampling logs (<--log-dir>/<chrom>.<path>.guess.real.log).
    Expects a line output by `./guess_n_and_cenhap.py`.
        Best guess: use <n> haplotypes & sample cenhap = <cenhap>.
    Also uses lines like 
        Selected haplotype <name> with score <score>
    to get the haplotypes sampled in order.
- Alignment stats (<--aln-dir>/<chrom>.<path>.stats.log).
    A stats file with lines
        [prefix]: identity [0-1] / correctness [%] / runtime [sec] / memory [GB]

== How I ran this ==

From within /private/home/fokamoto/centromere-haplotype-sampling-pipeline

./data_scripts/haploid_data.py \
    -c input_data -d input_data \
    -l /private/groups/patenlab/fokamoto/centrolign/graph/haploid \
    -a /private/groups/patenlab/fokamoto/centrolign/alignments/haploid > haploid_data.tsv
"""

import argparse # Command line argument parsing
import os # File system interaction
import sys # stderr
from typing import Dict, List, Tuple # Type hints

CenhapTable = Dict[str, str]
"""A {hap : cenhap} truth table."""
DistMatrix = Dict[str, Dict[str, float]]
"""A symmetric {hap1 : {hap2 : dist}} matrix."""

# Parts of file names (see file docstring)
CENHAP_SUFFIX = 'cenhap_predictions.tsv'
DIST_SUFFIX = 'r2_QC_v2_centrolign_pairwise_distance.csv'
GUESS_SUFFIX = 'guess.real.log'
STATS_SUFFIX = 'stats.log'

REFS = {'CHM13'      : 'CHM13', 
        'Neighbor'   : 'neighbor', 
        'Sampled'    : 'sampled', 
        'Native hap' : 'own_hap'}
"""Match up column titles to filenames"""
REALNESS = ['real', 'sim']
"""Possible levels of realness."""
ALIGNERS = ['minimap2', 'giraffe']
"""Names of aligners used."""

# All combinations, except that minimap2 can't align to the sampled graph
ALN_COMBOS = [(ref, realness, tool)
               for ref in REFS.keys()
               for realness in REALNESS
               for tool in ALIGNERS 
               if not (ref == 'Sampled' and tool == 'minimap2')]

def parse_args() -> argparse.Namespace:
    """Handle command-line argument parsing."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cenhap-dir', required=True,
                        help='Directory with TSVs of truth cenhap assignments')
    parser.add_argument('-d', '--dist-dir', required=True,
                        help='Directory with CSVs of pairwise distances')
    parser.add_argument('-l', '--log-dir', required=True,
                        help='Directory with haplotype sampling logs')
    parser.add_argument('-a', '--aln-dir', required=True,
                        help='Directory with read alignment stat files')
    return parser.parse_args()

def read_chrom_cenhap_table(cenhap_file: str) -> CenhapTable:
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

def read_all_cenhap_tables(cenhap_dir: str) -> Dict[str, CenhapTable]:
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

def read_chrom_distances(distances_file: str) -> DistMatrix:
    """Read a pairwise distance file.

    Three-column CSV with hap1,hap2,dist
    into symmetrical dictionary of
    {hap1 : {hap2 : dist}}.
    """
    dist_matrix = dict()
    with open(distances_file) as file:
        for line in file:
            parts = line.strip().split(',')
            hap1, hap2, dist = parts[0], parts[1], float(parts[2])

            # In case these haplotypes haven't been seen before
            if hap1 not in dist_matrix:
                dist_matrix[hap1] = dict()
            if hap2 not in dist_matrix:
                dist_matrix[hap2] = dict()
            
            dist_matrix[hap1][hap2] = dist
            dist_matrix[hap2][hap1] = dist
    return dist_matrix

def read_all_distances(dist_dir: str) -> Dict[str, DistMatrix]:
    """Reads all distance matrices in a directory.
    
    Finds all files matching
        <--dist-dir>/<chrom>_r2_QC_v2_centrolign_pairwise_distance.csv
    and reads them into a {chrom : matrix} dictionary
    """
    dist_matrices = dict()
    for item in os.listdir(dist_dir):
        if DIST_SUFFIX in item:
            chrom = item.split('_')[0]
            dist_matrices[chrom] = read_chrom_distances(
                os.path.join(dist_dir, item))
    return dist_matrices

def get_guesses(log_file: str) -> Tuple[List[str], List[str], int, str]:
    """Look up the sampled haplotypes & cenhap.
    
    Pulls specific sampled haplotypes via
        Selected haplotype <name> with score <score>
    And then n & cenhap via
        Best guess: use <n> haplotypes & sample cenhap = <cenhap>.

    Returns
    - A list of the haplotypes selected, in order
    - A list of their scores, in order
    - How many haplotypes from that were subsampled for the graph
    - What cenhap this haplotype was guessed as
    """

    sampled_haps = []
    scores = []
    n_haps = None
    guess_cenhap = None
    with open(log_file) as file:
        for line in file:
            parts = line.strip().split()
            if line.startswith('Selected haplotype'):
                hap_name = parts[2]
                sampled_haps.append(hap_name)
                scores.append(parts[5])
            elif line.startswith('Best guess'):
                n_haps = int(parts[3])
                guess_cenhap = parts[9]
            
    return sampled_haps, scores, n_haps, guess_cenhap

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

def write_data(cenhap_tables: Dict[str, CenhapTable],
               distance_matrices: Dict[str, DistMatrix],
               log_dir: str, aln_dir: str) -> None:
    """Write the output TSV (see file docstring)."""
    # Write header
    column_titles = ['Chromosome', 'Haplotype name',
                     'Truth cenhap', 'Guessed cenhap',
                     'Sampled haplotype names', 'Sampled haplotype scores',
                     'Minimum graph distance', 'Minimum sampled distance']
    for (ref, realness, tool) in ALN_COMBOS:
        aln_group = f'{ref} {tool} {realness}'
        column_titles += [f'{aln_group} identity', f'{aln_group} correctness',
                          f'{aln_group} runtime', f'{aln_group} memory']
    print('\t'.join(column_titles))

    for chrom in dist_matrices.keys():
        correct = 0
        total = 0
        for hap_name in dist_matrices[chrom].keys():
            if chrom in cenhap_tables and hap_name in cenhap_tables[chrom]:
                truth_cenhap = cenhap_tables[chrom][hap_name]
            else:
                truth_cenhap = None
            items_to_write = [chrom, hap_name, truth_cenhap]
            
            guess_file = os.path.join(
                log_dir, f'{chrom}.{hap_name}.{GUESS_SUFFIX}')
            if not os.path.exists(guess_file):
                # This haplotype didn't go through the pipeline; skip
                continue

            # Add guesses from logfile
            sampled_haps, scores, n_haps, guess_cenhap = get_guesses(guess_file)
            items_to_write += [guess_cenhap, ','.join(sampled_haps),
                               ','.join(scores)]

            # Look up distances
            dist_row = distance_matrices[chrom][hap_name]
            # Min distance to any haplotype
            true_closest_hap = min(dist_row, key=dist_row.get)
            # Min distance to a sampled haplotype
            closest_sampled_hap = min(sampled_haps[:n_haps], key=dist_row.get)
            items_to_write += [dist_row[true_closest_hap], 
                               dist_row[closest_sampled_hap]]

            # Dump in alignment stats as well
            aln_stats = get_aln_stats(os.path.join(
                aln_dir, f'{chrom}.{hap_name}.{STATS_SUFFIX}'))
            if len(ALN_COMBOS) != len(aln_stats):
                continue
            for (ref, realness, tool) in ALN_COMBOS:
                prefix = f'{chrom}.{hap_name}.{REFS[ref]}.{realness}.{tool}'
                items_to_write += aln_stats[prefix]

            if truth_cenhap == guess_cenhap and not truth_cenhap is None:
                correct += 1
            total += 1

            print('\t'.join(str(item) for item in items_to_write))
        print(f'{chrom} accuracy: {(100 * correct / total):.2f}%', 
              file=sys.stderr)

if __name__ == '__main__':
    args = parse_args()
    # Load cross-sample files first
    cenhap_tables = read_all_cenhap_tables(args.cenhap_dir)
    dist_matrices = read_all_distances(args.dist_dir)
    # Write by-sample data
    write_data(cenhap_tables, dist_matrices, args.log_dir, args.aln_dir)