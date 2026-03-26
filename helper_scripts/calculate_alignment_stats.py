#!/usr/bin/env python3
"""Calculate identity & correctness stats for reads.

== Inputs ==

- GFA file with haplotype paths (--gfa).
    Expects all paths as P lines; ignores all other lines.
    The node list is allowed to end with a comma.
    https://gfa-spec.github.io/GFA-spec/GFA1.html
- Truth read positions (--truth-tsv).
    Expects columns for read name and then comma-separated nodes.
    The node list is allowed to end with a comma.
    Also expects a header line but it will be skipped.
- Alignment metadata (--aln-tsv).
    Expects columns for read name, identity, and comma-separated nodes.
    The node list is allowed to end with a comma.
    Also expects a header line but it will be skipped.

== Correctness calculations ==

The typical definition of alignment correctness (i.e. position overlap)
is problematic for these "leave-one-out" alignments. If the "correct"
starting position is located in sequence private to that haplotype,
then the read would have nowhere "correct" to map to.

Thus, "correctness" is herein defined as % of possible node overlaps found.
That is, if the truth position covered 150 nodes, 50 of those nodes are on
private sequence, and the alignment overlaps 80 of the rest, then it has

    100% * 80 / (150 - 50) = 100% * 80 / 100 = 80% correctness

Orientation is ignored, and partial coverage is counted as full.
"""

import argparse # Command line argument parsing
from typing import Dict, Set # Type hints

def parse_args() -> argparse.Namespace:
    """Handle command-line argument parsing."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gfa', help='GFA file with path nodes')
    parser.add_argument('-p', '--path-name', help='Path reads are from')
    parser.add_argument('-t', '--truth-tsv', help='Truth nodes for reads')
    parser.add_argument('-a', '--aln-tsv', help='Aligned nodes for reads')
    parser.add_argument('-r', '--require-private', action='store_true',
                        help='Correctness calculations use private nodes')
    return parser.parse_args()

def parse_node_list(node_str: str) -> Set[int]:
    """Convert comma-separated list of node IDs to set of normalized IDs."""
    return {int(node.rstrip('+-')) for node 
            in node_str.strip().split(',') if node}

def find_private_nodes(gfa_file: str, path_name: str) -> Set[int]:
    """Use GFA file to find nodes private to a path."""
    on_this_path = set()
    on_other_paths = set()

    with open(gfa_file) as f:
        for line in f:
            if not line.startswith('P'):
                continue

            parts = line.rstrip().split()
            # Extract contig name from <sample ID>#<hap num>#<contig>
            cur_path_name = parts[1].split('#')[-1]
            nodes = parse_node_list(parts[2])

            if cur_path_name == path_name:
                on_this_path = nodes
            else:
                on_other_paths |= nodes
    
    # Subset to only private nodes
    return (on_this_path - on_other_paths)

def load_truth_nodes(read_tsv: str) -> Dict[str, Set[int]]:
    """Load the truth node positions for read set.
    
    Uses a TSV output by vg filter --tsv-out "name;nodes".
    Ignores orientations (see "Correctness calculations").
    Outputs {read name : {nodes}}
    """

    truth_nodes = dict()
    with open(read_tsv) as file:
        file.readline() # header
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                truth_nodes[parts[0]] = parse_node_list(parts[1])

    return truth_nodes

def calc_aln_stats(truth_nodes: Dict[str, Set[int]], private_nodes: Set[int],
                   aln_tsv_file: str, req_private: bool) -> None:
    """Calculate mean identity & correctness for an alignment set.
    
    Takes in an alignment TSV from vg filter --tsv-out "name;identity;nodes"
    and uses truth/private node lists to get correctness stats.
    Prints "mean identity = <identity> and mean correctness = <correctness>"

    req_private indicates whether these alignments should have used the
    native haplotype; if not, then no penalty is applied when nodes
    private to the haplotype path are missed. See "Correctness calculations".
    """

    identity = []
    correctness = []

    with open(aln_tsv_file) as file:
        file.readline() # header
        for line in file:
            parts = line.strip().split()
            identity.append(float(parts[1]))

            # Calculate correctness of this node list
            if len(parts) == 2 and parts[1] == '0':
                # Unmapped read has 0 correctness
                correctness.append(0)
            elif len(parts) == 3:
                aln_nodes = parse_node_list(parts[2])
                if parts[0] in truth_nodes:
                    truth = truth_nodes[parts[0]]

                    if not req_private:
                        # Don't require alignments to private nodes
                        truth -= private_nodes
                    if not truth:
                        # All truth nodes have been thrown out
                        continue
                    
                    overlap = truth & aln_nodes
                    correctness.append(100 * len(overlap) / len(truth))
            else:
                raise ValueError(f'No nodes in line from {aln_tsv_file}')

    print(f'mean identity = {sum(identity) / len(identity)} and '
          f'mean correctness = {sum(correctness) / len(correctness)}')

if __name__ == '__main__':
    args = parse_args()
    private_nodes = find_private_nodes(args.gfa, args.path_name)
    truth_nodes = load_truth_nodes(args.truth_tsv)
    calc_aln_stats(truth_nodes, private_nodes, args.aln_tsv, args.require_private)