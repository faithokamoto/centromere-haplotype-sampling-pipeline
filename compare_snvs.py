#!/usr/bin/env python3
"""Compate SNV calls with truth set and compute stats.

Input a VCF file of SNV calls and truth set CSV files.
- VCF file should be for a haploid sample and VCFv4.2 format
- Truth set CSVs uses columns 4, 6, and 7 from:
    ref_id,qry_id,var_type,ref_pos,qry_pos,ref_base,qry_base

One truth set CSV is the one that is actually used
and the other is relaxed to check for filtered out variants.

Prints out these stats:
- True Positives (TP): number of SNVs in both call set and truth set
- False Positives (FP): number of SNVs in call set but not in truth set
    - Also split this by whether they are in the relaxed truth set
- False Negatives (FN): number of SNVs in truth set but not in call set
    - Also split this by if I filtered it out or it was in my SV calls
- Precision: TP / (TP + FP)
- Recall: TP / (TP + FN)
"""

import argparse
from dataclasses import dataclass
from typing import List

@dataclass
class Variant:
    """A variant call, either from VCF or truth set."""
    ref_pos: int
    ref_allele: str
    alt_allele: str
    filtered: bool = False

    def is_snv(self) -> bool:
        """Does this variant represent an SNV?"""
        return len(self.ref_allele) == 1 and len(self.alt_allele) == 1
    
    def is_same_var(self, other) -> bool:
        """Are these two variants the same (regardless of type)?"""
        return (self.ref_pos == other.ref_pos and
                self.ref_allele == other.ref_allele and
                self.alt_allele == other.alt_allele)

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare SNV calls with truth set and compute stats."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="VCF file of SNV calls.",
    )
    parser.add_argument(
        "--truth-csv",
        required=True,
        help="Truth set CSV file.",
    )
    parser.add_argument(
        "--relaxed-truth-csv",
        required=True,
        help="Relaxed truth set CSV file for near-miss checking.",
    )
    return parser.parse_args()

def load_truth_set(csv_file: str) -> List[Variant]:
    """Load ordered truth set from CSV file."""
    truth_set = []
    with open(csv_file, 'r') as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.strip().split(',')
            var_type = parts[2]
            ref_pos = parts[3]
            ref_base = parts[5]
            qry_base = parts[6]

            # We don't care about truth-set SVs
            if var_type == "SNV":
                truth_set.append(Variant(int(ref_pos), ref_base, qry_base))
    return truth_set

def load_vcf(vcf_file: str) -> List[Variant]:
    """Load variant calls from VCF file."""
    call_set = []
    with open(vcf_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            # Convert to 0-based & get rid of dummy base
            ref_pos = int(parts[1]) - 2  
            ref_base = parts[3]
            alt_bases = parts[4].split(',')
            filtered = parts[6] != "PASS"

            # Handle multiple ALT alleles
            for alt_base in alt_bases:
                call_set.append(Variant(ref_pos, ref_base, alt_base, filtered))
    return call_set

def compute_stats(call_set: List[Variant], truth_set: List[Variant], 
                  relaxed_truth_set: List[Variant]) -> None:
    """Compute and print comparison stats.
    
    Take advantage of the fact that the Variants are in order;
    we can go through the lists and pop off elements as we find matches.
    """

    # The stats we want to track
    tp = 0
    fp = 0
    fn = 0
    fp_filtered = 0
    fn_filtered = 0
    fn_in_sv = 0

    # Pointers to current indexes
    truth_idx = 0
    call_idx = 0

    # Relaxed truth set is only used for looking up; {pos: Variant}
    relaxed_truth_dict = {var.ref_pos: var for var in relaxed_truth_set}

    while call_idx < len(call_set) and truth_idx < len(truth_set):
        call = call_set[call_idx]
        truth = truth_set[truth_idx]

        if not call.is_snv():
            # Not an SNV, so this doesn't count as a call
            # However, it may explain some FNs
            call_end = call.ref_pos + len(call.ref_allele) - 1

            if truth.ref_pos < call.ref_pos:
                # Any truth variants before this call are plain FNs
                fn += 1
                truth_idx += 1
            elif truth.ref_pos <= call_end:
                # Any truth variants covered by this call are FNs in SV
                fn_in_sv += 1
                fn += 1
                truth_idx += 1
            else:
                # Move to next call
                call_idx += 1
        else:
            # This is an SNV call
            if call.ref_pos < truth.ref_pos:
                # Call not in truth set, a FP
                # Check for near-miss in relaxed truth set
                if (call.ref_pos in relaxed_truth_dict
                    and relaxed_truth_dict[call.ref_pos].is_same_var(call)):
                    fp_filtered += 1
                    
                fp += 1
                call_idx += 1
            elif call.ref_pos > truth.ref_pos:
                # Truth variant not in calls
                fn += 1
                truth_idx += 1
            else:
                if truth.is_same_var(call):
                    # Match found; either TP or FN if filtered
                    if call.filtered:
                        fn_filtered += 1
                    else:
                        tp += 1
                else:
                    # Positions match but alleles don't; count as FN + FP
                    fn += 1
                    fp += 1
                call_idx += 1
                truth_idx += 1
    
    # Any remaining calls are FPs
    while call_idx < len(call_set):
        call = call_set[call_idx]
        if call.is_snv():
            # Check for being in relaxed truth set
            if (call.ref_pos in relaxed_truth_dict 
                and relaxed_truth_dict[call.ref_pos].is_same_var(call)):
                fp_filtered += 1
                
            fp += 1
        call_idx += 1

    # Any remaining truth variants are FNs
    while truth_idx < len(truth_set):
        fn += 1
        truth_idx += 1

    # Print stats
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    
    print(f"True Positives (TP): {tp}")
    print(f"False Positives (FP): {fp} (Filtered: {fp_filtered})")
    print(f"False Negatives (FN): {fn} (Filtered: {fn_filtered}, In SVs: {fn_in_sv})")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    # Also calculate recall with in-SV FNs removed
    recall_no_sv = tp / (tp + fn - fn_in_sv)
    print(f"Recall excluding FNs in SVs: {recall_no_sv:.4f}")

if __name__ == "__main__":
    args = parse_args()

    truth_set = load_truth_set(args.truth_csv)
    print(f"Loaded {len(truth_set)} SNVs from truth set.")
    relaxed_truth_set = load_truth_set(args.relaxed_truth_csv)
    print(f"Loaded {len(relaxed_truth_set)} SNVs from relaxed truth set.")
    call_set = load_vcf(args.vcf)
    print(f"Loaded {len(call_set)} variant calls from VCF.")

    compute_stats(call_set, truth_set, relaxed_truth_set)