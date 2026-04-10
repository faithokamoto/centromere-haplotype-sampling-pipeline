#!/usr/bin/env python3
"""Calculate depth on nodes private to each haplotype.

== Inputs ==

- GFA file with haplotype paths (--gfa).
    Expects all haplotypes as W lines; ignores all other lines.
    https://gfa-spec.github.io/GFA-spec/GFA1.html
- Alignment metadata (--aln-tsv).
    Expects columns for read name, identity, and comma-separated nodes.
    The node list is allowed to end with a comma.
    Also expects a header line but it will be skipped.
- Haplotype sampling log (--sampling-log)
    Expects a line output by `./guess_n_and_cenhap.py`.
        Best guess: use <n> haplotypes & sample cenhap = <cenhap>.
    Also uses lines like 
        Selected haplotype <name> with score <score>
    to get the haplotypes sampled in order.

== Private nodes ==

First we extract the set of haplotypes used in the personalized graph.
Then we look up their nodes in the GFA. A "private" node is one which
a given haplotype traverses but no other haplotype does within the
personalized graph. We then calculate average depth on these nodes.
"""

import argparse # Command line argument parsing
import re # Parsing walks (annoying format)
from typing import Dict, List, Set # Type hints

def parse_args() -> argparse.Namespace:
    """Handle command-line argument parsing."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gfa', help='GFA file with path nodes')
    parser.add_argument('-a', '--aln-tsv', help='Alignment TSV')
    parser.add_argument('-l', '--sampling-log', help='Haplotype sampling log')
    return parser.parse_args()

def parse_node_list(node_str: str) -> Set[int]:
    """Convert comma-separated list of node IDs to set of normalized IDs."""
    return {int(node.rstrip('+-')) for node 
            in node_str.strip().split(',') if node}

def parse_walk(walk: str) -> Set[int]:
    """Parse a W-line walk string like >12<34>56 into node IDs"""
    nodes = set()
    for match in re.finditer(r'[><](\d+)', walk):
        nodes.add(int(match.group(1)))
    return nodes

def read_nodes(gfa_file: str) -> Dict[str, Set[int]]:
    """Read nodes for each haplotype from W lines in a GFA"""
    hap_nodes = {}

    with open(gfa_file) as f:
        for line in f:
            if not line.startswith('W'):
                continue

            parts = line.rstrip().split()
            sample_id = parts[1]
            hap_index = parts[2]
            walk = parts[6]
            hap_nodes[f"{sample_id}.{hap_index}"] = parse_walk(walk)

    return hap_nodes

def find_private_nodes(nodes: Dict[str, Set[int]], hap_name: str) -> Set[int]:
    """Find which nodes are private to a particular haplotype."""
    seen_elsewhere = set()
    for other in nodes.keys():
        if other != hap_name:
            seen_elsewhere |= nodes[other]
    
    return (nodes[hap_name] - seen_elsewhere)

def get_guesses(log_file: str) -> List[str]:
    """Look up the sampled haplotypes.
    
    Pulls specific sampled haplotypes via
        Selected haplotype <name> with score <score>
    And grabs the exact number used via
        Best guess: use <n> haplotypes & sample cenhap = <cenhap>

    Returns the selected haplotypes in order
    """

    sampled_haps = []
    n_haps = None
    with open(log_file) as file:
        for line in file:
            parts = line.strip().split()
            if line.startswith('Selected haplotype'):
                sampled_haps.append(parts[2])
            elif line.startswith('Best guess'):
                n_haps = int(parts[3])
            
    return sampled_haps[:n_haps]

if __name__ == '__main__':
    args = parse_args()
    all_nodes = read_nodes(args.gfa)
    sampled_haps = get_guesses(args.sampling_log)
    all_nodes = {hap: nodes for hap, nodes in all_nodes.items()
                 if hap in set(sampled_haps)}
    
    # Precompute private nodes for top hap
    top_private = find_private_nodes(all_nodes, sampled_haps[0])
    top_coverage = 0

    with open(args.aln_tsv) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                # Unaligned reads have no node list; skip
                continue
            for node in parts[2].split(',')[:-1]:  # Ignore trailing comma
                node_id = int(node.rstrip('+-'))
                if node_id in top_private:
                    top_coverage += 1

    print(f'{sampled_haps[0]} has avg depth {top_coverage / len(top_private)}')