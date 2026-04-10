#!/usr/bin/env python3
"""For an input sample, guess *n* haps to sample & its cenhap.

Takes a haplotype sampling logfile and tries to guess
how many haplotypes should be sampled based on score.

Uses logfile lines of the form "Selected haplotype X with score Y".
Also takes a cenhap assignment file and tries to guess optimal cenhap.
"""

import argparse # Command-line argument parsing
from dataclasses import dataclass # Basic structs
from typing import Dict, List # Type hints

@dataclass
class SampledHaplotype:
    """Metadata about a haplotype from the sampler."""
    name: str
    score: float

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Guess n to sample & cenhap")
    parser.add_argument("--fall-threshold", '-f', type=int, default=500,
                        required=True, help="Limit for drop from top score")
    parser.add_argument("--cenhap-table", type=str, required=True,
                        help="Cenhap assignments")
    parser.add_argument("--dist-matrix", type=str, required=True,
                        help="Distance matrix")
    parser.add_argument("--ploidy", type=int, default=1, required=True,
                        help="Sample ploidy (1 or 2)")
    parser.add_argument("logfile", type=str, help="Logfile to read from")
    return parser.parse_args()

def parse_scores(logfile: str) -> List[SampledHaplotype]:
    """Extract scores from logfile.
    
    Reads 'Selected haplotype X with score Y' lines
    into list of SampledHaplotypes, assuming order.
    """
    sampled = []
    with open(logfile) as f:
        for line in f:
            # Extract X,Y from 'Selected haplotype X with score Y'
            if line.startswith('Selected haplotype'):
                parts = line.strip().split()
                sampled.append(SampledHaplotype(parts[2], float(parts[5])))
    return sampled

def read_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read the cenhap assignment table.
    
    Assuming the first line is a header,
    reads haplotype,cenhap lines into 
    {haplotype : cenhap} dictionary.
    """
    cenhap_table = {}
    with open(cenhap_file, 'r') as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            cenhap_table[parts[0]] = parts[1]
    return cenhap_table

def read_dist_matrix(matrix_file: str) -> Dict[str, Dict[str, float]]:
    """Read the distance matrix file.

    Three-column CSV with hap1,hap2,dist
    into symmetrical dictionary of
    {hap1 : {hap2 : dist}}.
    """
    dist_matrix = {}
    with open(matrix_file, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            hap1, hap2, dist = parts[0], parts[1], float(parts[2])

            if hap1 not in dist_matrix:
                dist_matrix[hap1] = {}
            if hap2 not in dist_matrix:
                dist_matrix[hap2] = {}
            
            dist_matrix[hap1][hap2] = dist
            dist_matrix[hap2][hap1] = dist
    return dist_matrix

def print_dist_info(dist_matrix: Dict[str, Dict[str, float]], 
                    sampled: List[SampledHaplotype], path_name: str) -> None:
    """Print distance context from matrix.

    For each sampled haplotype, reports its distance from truth.
    Also reports the top N closest haps, for N = |sampled|.
    """

    print("Sampled haplotypes:")

    for other_hap in sampled:
        other_name = other_hap.name
        print(f"{other_name} is dist {dist_matrix[path_name][other_name]}")

    print("Top haplotypes:")
    
    top_dist = sorted(dist_matrix[path_name], key=dist_matrix[path_name].get)
    for i in range(len(sampled)):
        top_name = top_dist[i]
        print(f"{top_name} is dist {dist_matrix[path_name][top_name]}")

def guess_optimal_n(sampled: List[SampledHaplotype], threshold: int) -> int:
    """Guess the optimal number of haplotypes to sample.

    If score drops more than threshold from the top,
    don't use the next sampled haplotype.
    """

    max_score = sampled[0].score
    
    for i in range(1, len(sampled)):
        cur_score = sampled[i].score

        # Check if we should stop here (= use prev i haps)
        if cur_score < max_score - threshold:
            # We've fallen too far
            return i
        
    # Guess we're using all of these
    return len(sampled)

def guess_cenhap(sampled: List[SampledHaplotype], cenhap_table: Dict[str, str],
                 n_to_use: int, ploidy: int) -> str:
    """Guess a sample's cenhap.
    
    For haploid samples, just use the top sample's cenhap.
    For diploid samples, keep going down the list until
    getting two different cenhaps, or report a homozygote.
    """

    if ploidy == 1:
        return cenhap_table.get(sampled[0].name)
    elif ploidy == 2:
        # Find cenhaps in use
        cenhaps_sampled = set()
        for hap in sampled[:n_to_use]:
            if hap.name in cenhap_table:
                cenhaps_sampled.add(cenhap_table.get(hap.name))
            if len(cenhaps_sampled) == 2:
                break
        
        # Construct genotype
        if not cenhaps_sampled:
            raise ValueError('No cenhaps sampled')
        elif len(cenhaps_sampled) == 1:
            cenhap = cenhaps_sampled.pop()
            return f'{cenhap} / {cenhap}'
        else:
            cenhap1 = cenhaps_sampled.pop()
            cenhap2 = cenhaps_sampled.pop()
            return f'{cenhap1} / {cenhap2}'
    
if __name__ == "__main__":
    args = parse_args()

    if args.ploidy != 1 and args.ploidy != 2:
        raise ValueError('--ploidy must be 1 or 2')

    scores = parse_scores(args.logfile)
    matrix = read_dist_matrix(args.dist_matrix)
    if args.ploidy == 1:
        hap_name_parts = args.logfile.split("/")[-1].split(".")
        hap_name = f'{hap_name_parts[1]}.{hap_name_parts[2]}'
        print_dist_info(matrix, scores, hap_name)
    optimal_n = guess_optimal_n(scores, args.fall_threshold)
    
    cenhap_table = read_cenhap_table(args.cenhap_table)
    cenhap = guess_cenhap(scores, cenhap_table, optimal_n, args.ploidy)
    
    print(f"Best guess: use {optimal_n} haplotypes & sample cenhap = {cenhap}")