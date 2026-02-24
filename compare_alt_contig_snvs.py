#!/usr/bin/env python3
"""Compate SNV calls with truth set and compute stats.

Input a VCF file of SNV calls and truth set CSV files,
plus an augmented reference segment map.
- VCF file should be for a haploid sample and VCFv4.2 format
- Truth set CSVs uses columns 4, 6, and 7 from:
    ref_id,qry_id,var_type,ref_pos,qry_pos,ref_base,qry_base
- Segment map TSV uses columns 1, 2, and 4 from (headerless):
    source_path,start_index,end_index,augref_path,ref_contig,ref_start,ref_end

Searches within a directory that holds a bunch of truth set CSVs.
Uses <ref>_<query>.snvs.500bp_95pct.csv

Prints out these stats:
- True Positives (TP): number of SNVs in both call set and truth set
- False Positives (FP): number of SNVs in call set but not in truth set
- Precision: TP / (TP + FP)
"""

import argparse # Command line argument parsing
from dataclasses import dataclass # Simple C++ like struct
from typing import Dict, Set # Type hints

TSV_SUFFIX = 'snvs.10bp_95pct.csv'
CHM13_NAME = 'CHM13.0'
WIGGLE_ROOM = 5

@dataclass
class SNV:
    """An SNV call, either from VCF or truth set."""
    ref_pos: int
    ref_allele: str
    alt_allele: str

    def __hash__(self):
        return hash((self.ref_pos, self.ref_allele, self.alt_allele))

@dataclass
class AltContig:
    """An alternative contig path from a non-reference haplotype."""
    source_path: str
    # Indexes of the start of this segment within the source path
    start_index: int

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Compare SNV calls with truth set and compute stats.'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        help='VCF file of SNV calls.',
    )
    parser.add_argument(
        '-s',
        '--segment-map',
        help='Segment map TSV file for augref contig names.'
    )
    parser.add_argument(
        '-d',
        '--truth-dir',
        help='Truth set CSV directory.',
    )
    parser.add_argument(
        '-p',
        '--path-name',
        help='Name of path that calls are for.',
    )
    return parser.parse_args()

def load_truth_set(csv_file: str) -> Set[SNV]:
    """Load ordered truth set from CSV file.
    
    Convert a CSV with header row
        ref_id,qry_id,var_type,ref_pos,qry_pos,ref_base,qry_base
    to a set of SNVs which store (ref_pos, ref_base, qry_base)

    Any SVs are filtered out.
    """

    truth_set = set()
    with open(csv_file) as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.strip().split(',')
            ref_pos = int(parts[3])
            ref_base = parts[5]
            qry_base = parts[6]

            # Skip SVs
            if len(ref_base) != 1 or len(qry_base) != 1 or qry_base == '*':
                continue

            truth_set.add(SNV(ref_pos, ref_base, qry_base))
    return truth_set

def load_vcf(vcf_file: str) -> Dict[str, Set[SNV]]:
    """Load SNV calls from VCF file.
    
    Convert a v4.2 VCF for a single sample with header
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
    to a dictionary of SNVs with
        {contig : {(ref_pos, ref_base, alt_base)}}
    
    As these are used as variant calls, the implicit assumption
    is that all variants in the file are called as alts.
    """

    call_set = dict()
    with open(vcf_file) as f:
        if f.readline().strip() != '##fileformat=VCFv4.2':
            raise ValueError(f'Wrong VCF version: not 4.2 but\n{line}')
        
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            ref_contig = parts[0]
            # Convert to 0-based & get rid of dummy base
            ref_pos = int(parts[1]) - 2
            ref_base = parts[3]
            alt_base = parts[4]
            filtered = parts[6] != 'PASS'

            # Skip things that the VCF filtered out
            if filtered:
                continue
            # Skip anything that has more than one alt allele
            if ',' in alt_base:
                continue
            # Skip SVs
            if len(ref_base) != 1 or len(alt_base) != 1 or alt_base == '*':
                continue

            if not ref_contig in call_set:
                call_set[ref_contig] = set()
            call_set[ref_contig].add(SNV(ref_pos, ref_base, alt_base))
    return call_set

def load_segment_map(tsv_file: str) -> Dict[str, AltContig]:
    """Load alternative contigs from an augref segment map.
    
    Convert a TSV without header row
        source_path,start_index,end_index,augref_path,ref_contig,ref_start,ref_end
    to a dictionary of {new_name : {(source_path, start_index)}}
    """

    segments = dict()
    with open(tsv_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            source_path = parts[0]
            start_index = int(parts[1])
            augref_name = parts[3]

            segments[augref_name] = AltContig(source_path, start_index)
    return segments

def report_stats(tp: int, fp: int) -> None:
    """Print TP, FP, and precision stats"""
    precision = tp / (tp + fp)
    print(f'True Positives (TP): {tp}')
    print(f'False Positives (FP): {fp}')
    print(f'Precision: {precision:.4f}')
    print()

def compute_stats(segment_map: Dict[str, AltContig],
                  call_set: Dict[str, Set[SNV]], 
                  truth_dir: str, path_name: str) -> None:
    """Compute and print comparison stats.

    For each ref contig, pull the relevant truth sets and
    calculate stats. Also aggergate stats across contigs.
    """

    # The stats we want to track
    tp = 0
    fp = 0
    alt_tp = 0
    alt_fp = 0

    for cur_contig, cur_calls in call_set.items():
        # Look up the correct reference contig name
        if cur_contig == CHM13_NAME:
            ref_name = CHM13_NAME
            offset = 0
        else:
            ref_contig = segment_map[cur_contig]
            # Extract haplotype name from path name
            ref_name = ref_contig.source_path.split('#')[-1]
            offset = ref_contig.start_index
            print(f'translated {cur_contig} to {ref_name}')
        
        truth_file = f'{truth_dir}/{ref_name}_{path_name}.{TSV_SUFFIX}'
        
        try:
            truth_set = load_truth_set(truth_file)
            pos_indexed_truth = {snv.ref_pos: snv for snv in truth_set}
            
            cur_tp = 0
            cur_fp = 0

            for call in cur_calls:
                call.ref_pos += offset
                if call in truth_set:
                    cur_tp += 1
                else:
                    cur_fp += 1
                    printed_call = False
                    for pos in range(call.ref_pos - WIGGLE_ROOM, 
                                     call.ref_pos + WIGGLE_ROOM + 1):
                        if pos in pos_indexed_truth:
                            if not printed_call:
                                print(f'SNVs nearby {call}:')
                                printed_call = True
                            print('\t' + str(pos_indexed_truth[pos]))

            print(f'Calls against {cur_contig}:')
            report_stats(cur_tp, cur_fp)

            tp += cur_tp
            fp += cur_fp
            if ref_name != CHM13_NAME:
                alt_tp += cur_tp
                alt_fp += cur_fp

        except FileNotFoundError:
            print(f'{len(cur_calls)} calls against {ref_name}; no truth CSV')
            print()

    print('Overall:')
    report_stats(tp, fp)
    print('On alt contigs:')
    report_stats(alt_tp, alt_fp)

if __name__ == '__main__':
    args = parse_args()
    call_set = load_vcf(args.vcf)
    segment_map = load_segment_map(args.segment_map)
    compute_stats(segment_map, call_set, args.truth_dir, args.path_name)