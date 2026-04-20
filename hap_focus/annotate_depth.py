#!/usr/bin/env python3
"""Add DP tags to a GFA once given depth per position."""

import argparse # Command-line argument parsing
from collections import defaultdict # Easier counting

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Annotate GFA nodes with average depth (DP tag).')
    parser.add_argument('-g', '--gfa', required=True, help='Input GFA file')
    parser.add_argument('-d', '--depths', required=True, help='Depths TSV')
    return parser.parse_args()

def load_depths(depths_file: str) -> dict:
    """Calculate average depth on each node"""
    depth_sums = defaultdict(float)
    depth_counts = defaultdict(int)

    with open(depths_file) as file:
        file.readline() # Skip header line
        for line in file:
            parts = line.strip().split('\t')
            node_id = parts[1]
            depth_sums[node_id] += int(parts[3])
            depth_counts[node_id] += 1

    return {node : depth_sums[node] / depth_counts[node] 
            for node in depth_sums.keys()}

def annotate_gfa(input_gfa: str, depth_avg: dict) -> None:
    """Output an GFA with DP tags added on."""
    with open(input_gfa) as file:
        for line in file:
            if line.startswith('S\t'):
                fields = line.strip().split('\t')
                node_id = fields[1]

                if node_id in depth_avg:
                    fields.append(f'DP:f:{depth_avg[node_id]:.3f}')

                print('\t'.join(fields))
            else:
                print(line.strip())

if __name__ == '__main__':
    args = parse_args()
    depth_avg = load_depths(args.depths)
    annotate_gfa(args.gfa, depth_avg)