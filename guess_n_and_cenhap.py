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
                        help="Hard limit for score drop from top score")
    parser.add_argument("--plateau-threshold", '-p', type=int, default=50,
                        help="Required score drop between haplotypes")
    parser.add_argument("--cenhap-table", type=str, help="Cenhap assignments")
    parser.add_argument("--ploidy", type=int, help="Sample ploidy (1/2)")
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
            parts = line.strip().split(',')
            haplo, cenhap = parts[0], parts[1]
            # Handle missing cenhap assignments
            if cenhap != '':
                cenhap_table[haplo] = cenhap
    return cenhap_table

def guess_optimal_n(sampled: List[SampledHaplotype], fall_threshold: int,
                    plateau_threshold: int) -> int:
    """Guess the optimal number of haplotypes to sample.

    If score drops more than fall_threshold from the top,
    or drops less than plateau_threshold between haps,
    don't use the next sampled haplotype.
    """

    max_score = sampled[0].score
    prev_score = sampled[0].score

    left_first_plateau = False
    
    for i in range(1, len(sampled)):
        cur_score = sampled[i].score

        if cur_score <= prev_score - plateau_threshold:
            left_first_plateau = True

        # Check if we should stop here (= use prev i haps)
        if cur_score < max_score - fall_threshold:
            # We've fallen too far
            return i
        elif left_first_plateau and cur_score > prev_score - plateau_threshold:
            # We're not falling fast enough (hit plateau)
            return i
        
        prev_score = cur_score
        
    # Guess we're using all of these
    return len(sampled)

def guess_cenhap(sampled: List[SampledHaplotype], cenhap_table: Dict[str, str], 
                 ploidy: int, n_to_use: int) -> str:
    """Guess the cenhap for the sampled haplotypes.
    
    Look up cenhap assignments for top n sampled haplotypes
    and use that to guess cenhap via weighted 'voting'
    (last sampled gets 1, next up gets 2, etc.)
    """
    # Calculate weighted score for each cenhap type
    cenhap_scores = {}
    for i in range(n_to_use):
        cenhap = cenhap_table.get(sampled[i].name, None)
        if cenhap is not None:
            # Lower i get bigger votes
            vote = n_to_use - i
            cenhap_scores[cenhap] = cenhap_scores.get(cenhap, 0) + vote
    # Guess the cenhap with the highest score
    if not cenhap_scores:
        return 'Unknown'
    else:
        if ploidy == 2:
            # For diploid, return both cenhaps separated by '/'
            sorted_cenhaps = sorted(cenhap_scores.items(), 
                                    key=lambda x: x[1], reverse=True)
            if len(sorted_cenhaps) >= 2:
                return f"{sorted_cenhaps[0][0]}/{sorted_cenhaps[1][0]}"
            else:
                return f"{sorted_cenhaps[0][0]}/{sorted_cenhaps[0][0]}"
        else:
            return max(cenhap_scores.items(), key=lambda x: x[1])[0]

if __name__ == "__main__":
    args = parse_args()

    scores = parse_scores(args.logfile)
    optimal_n = guess_optimal_n(scores, args.fall_threshold,
                                args.plateau_threshold)
    
    cenhap_table = read_cenhap_table(args.cenhap_table)
    cenhap = guess_cenhap(scores, cenhap_table, args.ploidy, optimal_n)
    
    print(f"Best guess: use {optimal_n} haplotypes & sample cenhap = {cenhap}")