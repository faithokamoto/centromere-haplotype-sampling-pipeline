#!/usr/bin/env python3
"""Guess optimal *n* for how many haplotypes to sample.

Input a haplotype name and this script will automatically find the
relevant file locations. Output will be the value of *n* to use,
or 0 if the haplotype is deemed "hopeless".
"""

import argparse
import csv
import os

from typing import Dict, Set

SAMPLE_TABLE = 'input_data/test_samples.txt'
"""Four-column CSV of (path name, version, haplotype name, optimal N)."""
TSV_DIR = "/private/groups/patenlab/fokamoto/centrolign/alignments/leave_one_out"
"""Where alignments to haplotype-sampled graphs are stored."""

MAX_HOPELESS_FRACTION = 0.5
"""Max fraction of reads aligned with < HOPELESS_IDENTITY_THRESHOLD identity."""
HOPELESS_IDENTITY_THRESHOLD = 0.99
"""Identity threshold below which a read is considered poorly aligned."""
JUMP_IDENTITY_THRESHOLD = 0.0025
"""Average identity change threshold to consider a jump significant."""

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Guess optimal number of haplotypes to sample")
    parser.add_argument("--hap1-reads", type=str, required=False, default=None,
                        help="FASTQ file of reads from haplotype 1")
    parser.add_argument("--hap2-reads", type=str, required=False, default=None,
                        help="FASTQ file of reads from haplotype 2")
    parser.add_argument("name", type=str,
                        help="Name of haplotype/sample (e.g., 'HG02622_mat')")
    return parser.parse_args()

def collect_tsvs(name: str) -> Dict[int, str]:
    """Collect relevant TSV files for a given name.

    Parameters
    ----------
    name : str
        The haplotype name (e.g., 'HG02622_mat').

    Returns
    -------
    tsv_files : Dict[int, str]
        A dictionary mapping number of sampled haplotypes to TSV file paths.
    """

    tsv_files = dict()
    for _, __, files in os.walk(TSV_DIR):
        for file in files:
            if file.startswith(f"real_{name}") and file.endswith(".tsv"):
                # Parse real_<name>.<num_hap>haps.tsv
                parts = file.split(".")
                if len(parts) != 3:
                    continue
                num_hap_str = parts[1][0]
                num_hap = int(num_hap_str)
                tsv_files[num_hap] = os.path.join(TSV_DIR, file)

    return tsv_files

def read_identity_tsv(tsv_file: str, read_names: Set[str]) -> Dict[str, float]:
    """Read read identity from a TSV file.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file.
    read_names : Set[str]
        A set of read names to filter to.

    Returns
    -------
    identities : Dict[str, float]
        A dictionary mapping read names to their identity values.
    """
    
    identities = dict()
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            read_name = row['#name']
            if read_names is None or read_name in read_names:
                identity = float(row['identity'])
                identities[read_name] = identity
    return identities

def is_file_hopeless(identities: Dict[str, float]) -> bool:
    """Decide if a given TSV file indicates a hopeless case.
    
    Hopeless case means that < MAX_HOPELESS_FRACTION of reads
    aligned with >= HOPELESS_IDENTITY_THRESHOLD identity.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file.

    Returns
    -------
    is_hopeless : bool
        True if the case is hopeless, False otherwise.
    """

    total_reads = 0
    hopeless_reads = 0
    for identity in identities.values():
        total_reads += 1
        if identity < HOPELESS_IDENTITY_THRESHOLD:
            hopeless_reads += 1
    if total_reads == 0:
        return True
    hopeless_fraction = hopeless_reads / total_reads
    return hopeless_fraction > MAX_HOPELESS_FRACTION

def calculate_average_identity(identities: Dict[str, float]) -> float:
    """Calculate average identity for a set of alignments.

    Parameters
    ----------
    identities : Dict[str, float]
        A dictionary mapping read names to their identity values

    Returns
    -------
    average_identity : float
        The average identity value.
    """

    if len(identities) == 0:
        return 0.0
    total_identity = sum(identities.values())
    return total_identity / len(identities)

def guess_optimal_n(tsv_files: Dict[int, str], 
                    read_names: Set[str] = None) -> int:
    """Guess the optimal number of haplotypes to sample.

    Parameters
    ----------
    tsv_files : Dict[int, str]
        A dictionary mapping number of sampled haplotypes to TSV file paths.
    read_names : Set[str]
        A set of read names to filter to.

    Returns
    -------
    optimal_n : int
        The guessed optimal number of haplotypes to sample.
    """

    found_non_hopeless = False
    n_with_jump = None
    # Average identity for the alignment set which had a jump
    jump_avg = 0.0
    for n in sorted(tsv_files.keys()):
        current_identities = read_identity_tsv(tsv_files[n], read_names)
        cur_avg = calculate_average_identity(current_identities)
        print(f"n={n}: avg identity={cur_avg:.6f}")
        # Don't bother checking out this n further if hopeless
        if is_file_hopeless(current_identities):
            continue
        
        if not found_non_hopeless:
            found_non_hopeless = True
            # Going from hopeless to non-hopeless is a jump
            n_with_jump = n
            jump_avg = calculate_average_identity(current_identities)
        else:
            cur_avg = calculate_average_identity(current_identities)
            if cur_avg - jump_avg > JUMP_IDENTITY_THRESHOLD:
                n_with_jump = n
                jump_avg = cur_avg
                print(f"  Significant jump detected at n={n}!")
    
    return n_with_jump if found_non_hopeless else 0

def read_fastq_read_names(fastq_file: str) -> Set[str]:
    """Read read names from a FASTQ file."""

    read_names = set()
    with open(fastq_file, 'r') as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line_num % 4 == 1:  # Header line
                read_name = line.strip().split()[0][1:]  # Remove '@'
                read_names.add(read_name)
    return read_names

if __name__ == "__main__":
    args = parse_args()

    if (args.hap1_reads is not None) != (args.hap2_reads is not None):
        print("Error: Must provide both --hap1-reads and --hap2-reads or neither.")
        exit(1)

    tsv_files = collect_tsvs(args.name)
    if len(tsv_files) == 0:
        print(f"No TSV files found for name: {args.name}")
        exit(1)
    
    # Guess overall optimal number
    optimal_n = guess_optimal_n(tsv_files)
    if optimal_n == 0:
        print(f"{args.name} is hopeless")
    else:
        print(f"Sample {optimal_n} haplotypes for {args.name}")
    
    if args.hap1_reads is not None and args.hap2_reads is not None:
        print("Guessing n separately for each haplotype...")
        # Both haplotypes' reads provided; guess separately
        hap1_read_names = read_fastq_read_names(args.hap1_reads)
        hap2_read_names = read_fastq_read_names(args.hap2_reads)

        optimal_n_hap1 = guess_optimal_n(tsv_files, hap1_read_names)
        print("----")
        optimal_n_hap2 = guess_optimal_n(tsv_files, hap2_read_names)
        
        if optimal_n_hap1 == 0:
            print(f"{args.name} haplotype 1 is hopeless")
        else:
            print(f"Sample {optimal_n_hap1} haplotypes for haplotype 1")
        if optimal_n_hap2 == 0:
            print(f"{args.name} haplotype 2 is hopeless")
        else:
            print(f"Sample {optimal_n_hap2} haplotypes for haplotype 2")