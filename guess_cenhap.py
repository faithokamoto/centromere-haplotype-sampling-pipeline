#!/usr/bin/env python3
"""Guess sample cenhap.

Input a haplotype name and this script will attempt to guess
what cenhap it should be, based on the haplotyeps sampled at
the optimal *n*. If it is hopeless, it will say so.
"""

import argparse

from typing import Dict, List, Tuple

"""Four-column CSV of (path name, version, haplotype name, optimal N)."""
CENHAP_TABLE = '/private/groups/migalab/juklucas/centrolign/notes/correct_cenhaps_chr12/chr12_cenhap_assignments_final.csv'

SampledHaplo = List[Tuple[str, float]]
"""List of (haplotype, score) tuples."""

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Guess cenhap(s) of input")
    parser.add_argument("--ploidy", type=int, required=True, 
                        help="Ploidy of the sample")
    parser.add_argument("--logfile", type=str, required=True, 
                        help="Log file to read")
    parser.add_argument("--force-n", type=int, default=None,
                        help="Force target n instead of using log file")
    parser.add_argument("--max-n", type=int, default=None,
                        help="Maximum n to consider (only if forcing n)")
    return parser.parse_args()

def get_path_name(name: str, ploidy: int) -> set:
    """Get paths from a haplotype or sample name."""
    if '.' in name:
        # Already a path name
        return {name}
    
    if ploidy == 2:
        return {f"{name}.1", f"{name}.2"}
    elif ploidy == 1:
        if name.endswith('_pat') or name.endswith('_hap1'):
            return {name.split('_')[0] + '.1'}
        elif name.endswith('_mat') or name.endswith('_hap2'):
            return {name.split('_')[0] + '.2'}
    
    raise ValueError("Ploidy must be 1 or 2.")

def read_optimal_n(log_file: str, force_n: int = None) -> Tuple[str, int]:
    """Read name and optimal *n* from log file."""
    with open(log_file, 'r') as file:
        # Processing sample: <name>
        name = file.readline().strip().split()[-1]
        print("Processing sample:", name)
        if force_n is not None:
            return name, force_n
        else:
            # Read until we find the optimal n line
            for line in file:
                if f"haplotypes for {name}" in line:
                    # Sample <n> haplotypes for <name>
                    return name, int(line.strip().split()[1])
                if line.strip() == f"{name} is hopeless":
                    return name, 0
    raise ValueError(f"{name} not found in log file.")

def read_sampled_haplotypes(log_file: str, name: str, n_hap: int) -> SampledHaplo:
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
            if 'Processing sample' in line and name in line:
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
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split(',')
            haplo, cenhap = parts[0], parts[1]
            cenhap_table[haplo] = cenhap
    return cenhap_table

def guess_cenhap(sampled: SampledHaplo, cenhap_table: Dict[str, str], 
                 ploidy: int, n_to_use: int) -> str:
    """Guess the cenhap for the sampled haplotypes."""
    # Normalize scores so minimum is 1
    min_score = min(score for _, score in sampled)
    normalized = [(haplo, score - min_score + 1) for haplo, score in sampled]
    # Calculate weighted score for each cenhap type
    cenhap_scores = {}
    cenhaps_used = 0
    for haplo, score in normalized:
        cenhap = cenhap_table.get(haplo, None)
        if cenhap is not None:
            cenhaps_used += 1
            print("\tHaplotype:", haplo, "Cenhap:", 
                  cenhap, "Score:", score + min_score - 1)
            cenhap_scores[cenhap] = cenhap_scores.get(cenhap, 0) + score
            # Have we used enough haplotypes?
            if cenhaps_used == n_to_use:
                break
    # Guess the cenhap with the highest score
    if not cenhap_scores:
        return 'Unknown'
    else:
        if ploidy == 2:
            # For diploid, return both cenhaps separated by '/'
            sorted_cenhaps = sorted(cenhap_scores.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_cenhaps) >= 2:
                return f"{sorted_cenhaps[0][0]} {sorted_cenhaps[1][0]}"
            else:
                return f"{sorted_cenhaps[0][0]} {sorted_cenhaps[0][0]}"
        else:
            return max(cenhap_scores.items(), key=lambda x: x[1])[0]

if __name__ == '__main__':
    args = parse_args()

    if args.ploidy not in (1, 2):
        raise ValueError("Ploidy must be 1 or 2.")
    if args.force_n is not None and args.max_n is not None:
        if args.force_n > args.max_n:
            raise ValueError("force-n cannot be greater than max-n.")

    name, optimal_n = read_optimal_n(args.logfile, args.force_n)
    if optimal_n == 0:
        print(f"{name} has optimal n=0; cannot guess cenhap.")
        exit(0)
    n_to_read = optimal_n if args.max_n is None else args.max_n
    sampled_haplotypes = read_sampled_haplotypes(args.logfile, name, n_to_read)
    cenhap_table = read_cenhap_table(CENHAP_TABLE)

    guessed_cenhap = guess_cenhap(sampled_haplotypes, cenhap_table,
                                  args.ploidy, optimal_n)
    correct_cenhap = ''
    for path in get_path_name(name, args.ploidy):
        correct_cenhap += cenhap_table.get(path, 'Unknown') + ' '
    print(f"Guessed cenhap for {name} (optimal n={optimal_n}): "\
          f"{guessed_cenhap} (correct: {correct_cenhap.strip()})")