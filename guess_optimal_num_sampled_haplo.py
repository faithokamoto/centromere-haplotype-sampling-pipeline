#!/usr/bin/env python3
"""Guess optimal *n* for how many haplotypes to sample.

Input a haplotype name and this script will automatically find the
relevant file locations. Output will be the value of *n* to use,
or 0 if the haplotype is deemed "hopeless".
"""

import argparse
import csv
import os

from typing import Dict

BED_DIR = "/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files"
"""Mira's BED files, the names of which hint at path-to-haplotype mapping."""
TSV_DIR = "/private/groups/patenlab/fokamoto/centrolign/alignments/leave_one_out"
"""Where alignments to haplotype-sampled graphs are stored."""

MAX_HOPELESS_FRACTION = 0.5
"""Max fraction of reads aligned with < HOPELESS_IDENTITY_THRESHOLD identity."""
HOPELESS_IDENTITY_THRESHOLD = 0.99
"""Identity threshold below which a read is considered poorly aligned."""
JUMP_IDENTITY_THRESHOLD = 0.001
"""Average identity change threshold to consider a jump significant."""

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Guess optimal number of haplotypes to sample")
    parser.add_argument("haplotype", type=str,
                        help="Name of the haplotype (e.g., 'HG02622_mat')")
    return parser.parse_args()

def convert_path_to_haplotype(path_name: str) -> str:
    """Convert a path name to a haplotype name.

    In an attempt to be a little more user-friendly,
    this function automatically translates a path name
    (e.g. as you would find in the GFA) to the
    corresponding haplotype name used in file names.

    Parameters
    ----------
    path_name : str
        The path name (e.g., 'HG02622.2').

    Returns
    -------
    haplotype_name : str
        The haplotype name (e.g., 'HG02622_mat').
    """

    # Parse sections of <sample ID>.<haplotype number>
    if not "." in path_name:
        raise ValueError(f"Invalid path name: {path_name}")
    sample_id, hap_num = path_name.split(".")

    # If this were a pat/mat trio sample, what would the haplotype be?
    if hap_num == "1":
        name_if_trio = f"{sample_id}_pat"
    elif hap_num == "2":
        name_if_trio = f"{sample_id}_mat"
    else:
        raise ValueError(f"Invalid haplotype number in path name: {path_name}")

    # Check if the corresponding BED file exists
    for _, __, files in os.walk(BED_DIR):
        for file in files:
            if file.startswith(name_if_trio) and file.endswith(".bed"):
                return name_if_trio
    
    # Otherwise, this is a non-trio hap1/hap2 sample
    return f"{sample_id}_hap{hap_num}"

def collect_tsvs(haplotype: str) -> Dict[int, str]:
    """Collect relevant TSV files for a given path name.

    Parameters
    ----------
    haplotype : str
        The haplotype name (e.g., 'HG02622_mat').

    Returns
    -------
    tsv_files : Dict[int, str]
        A dictionary mapping number of sampled haplotypes to TSV file paths.
    """

    tsv_files = dict()
    for _, __, files in os.walk(TSV_DIR):
        for file in files:
            if file.startswith(f"real_{haplotype}") and file.endswith(".tsv"):
                # Parse real_<haplotype>.<num_hap>haps.tsv
                parts = file.split(".")
                if len(parts) != 3:
                    continue
                num_hap_str = parts[1][0]
                num_hap = int(num_hap_str)
                tsv_files[num_hap] = os.path.join(TSV_DIR, file)

    return tsv_files

def read_identity_tsv(tsv_file: str) -> Dict[str, float]:
    """Read read identity from a TSV file.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file.

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

def calculate_average_identity_change(identities1: Dict[str, float],
                                      identities2: Dict[str, float]) -> float:
    """Calculate average identity change between two sets of identities.

    Assumes that all reads are in both dictionaries,
    as these are supposed to be alignments of the
    same reads to different graphs.

    Parameters
    ----------
    identities1 : Dict[str, float]
        A dictionary mapping read names to their identity values (first set).
    identities2 : Dict[str, float]
        A dictionary mapping read names to their identity values (second set).

    Returns
    -------
    avg_change : float
        The average change in identity between the two sets.
    """

    total_change = 0.0
    count = 0
    for read_name in identities1.keys():
        change = identities2[read_name] - identities1[read_name]
        total_change += change
        count += 1
    if count == 0:
        return 0.0
    return total_change / count

def guess_optimal_n(tsv_files: Dict[int, str]) -> int:
    """Guess the optimal number of haplotypes to sample.

    Parameters
    ----------
    tsv_files : Dict[int, str]
        A dictionary mapping number of sampled haplotypes to TSV file paths.

    Returns
    -------
    optimal_n : int
        The guessed optimal number of haplotypes to sample.
    """

    found_non_hopeless = False
    n_with_jump = None
    previous_identities = None
    for n in sorted(tsv_files.keys()):
        current_identities = read_identity_tsv(tsv_files[n])
        # Don't bother checking out this n further if hopeless
        if is_file_hopeless(current_identities):
            continue
        
        if not found_non_hopeless:
            found_non_hopeless = True
            # Going from hopeless to non-hopeless is a jump
            n_with_jump = n
        else:
            avg_change = calculate_average_identity_change(
                previous_identities, current_identities)
            print(f"Avg identity delta from n={n-1} to n={n}: {avg_change:.6f}")
            if avg_change > JUMP_IDENTITY_THRESHOLD:
                n_with_jump = n
                print("JUMP!!")
        
        previous_identities = current_identities
    
    return n_with_jump if found_non_hopeless else 0

if __name__ == "__main__":
    args = parse_args()

    tsv_files = collect_tsvs(args.haplotype)
    if len(tsv_files) == 0:
        raise ValueError(f"No TSV files found for haplotype: {args.haplotype}")

    optimal_n = guess_optimal_n(tsv_files)
    if optimal_n == 0:
        print(f"Haplotype {args.haplotype} is hopeless")
    else:
        print(f"Sample {optimal_n} haplotypes for {args.haplotype}")