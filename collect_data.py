#!/usr/bin/env python3
"""Collect data about cenhap alignment attempts.

== Output ==

For each haplotype run, collect
- Path name
- Truth cenhap
- Guessed cenhap
- Num haplotypes sampled
- Distance to closest hap in graph
- Distance to closest sampled hap (within the ones we chose to align to)
- Identity & correctness for each alignment run:
    - Minimap2 CHM13
    - Giraffe CHM13
    - Minimap2 neighbor
    - Giraffe neighbor
    - Giraffe sampled
    - Minimap2 native hap
    - Giraffe native hap

These are output as a TSV to standard output.

== Inputs ==

General inputs, processed before specific haplotypes:
- Truth cenhap assignment TSV (--cenhap-table).
    Expects columns for path name and then cenhap, in that order.
    Also expects a header line but it will be skipped.
- Pairwise haplotype distances (--distances).
    Expects columns for hap1, hap2, and then float distance (0-1).
    Does not expect a header line, or even symmetric lines.
- GFA file with haplotype paths (--gfa).
    Expects all paths as P lines; ignores all other lines.
    The node list is allowed to end with a comma.
    https://gfa-spec.github.io/GFA-spec/GFA1.html

Sample-specific inputs:
- Haplotype sampling logs (<--log-dir>/<path>.guess.log).
    Expects a line output by `./guess_n_and_cenhap.py`.
    "Best guess: use <n> haplotypes & sample cenhap = <cenhap>".
    Also uses lines like "Selected haplotype <name> with score <score>"
    to get the haplotypes sampled in order.
- Truth read positions (<--reads-dir>/real_<path>.chr12.hifi.tsv).
    Expects columns for read name and then comma-separated nodes.
    The node list is allowed to end with a comma.
    Also expects a header line but it will be skipped.
- Alignment metadata (<--aln-dir>/<path>.<suffix>), with ALN_SUFFIXES suffixes.
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

The only variation to this formula occurs for "native hap" alignments,
where as the haplotype was not removed, the private nodes are also required.
"""

import argparse # Command line argument parsing
import os # File system interaction
from typing import Dict, List, Tuple, Set # Type hints

ALN_SUFFIXES = {
         'CHM13 minimap2' : 'chm13.real.minimap2.tsv',
          'CHM13 giraffe' : 'chm13.real.giraffe.tsv',
      'Neighbor minimap2' : 'neighbor.real.minimap2.tsv',
       'Neighbor giraffe' : 'neighbor.real.giraffe.tsv',
        'Sampled giraffe' : 'sampled.real.giraffe.tsv',
    'Native hap minimap2' : 'own_hap.real.minimap2.tsv',
     'Native hap giraffe' : 'own_hap.real.giraffe.tsv',
}
"""Alignment runs matched to their TSV suffixes."""

def parse_args() -> argparse.Namespace:
    """Handle command-line argument parsing."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cenhap-table',
                        help='TSV of truth cenhap assignments')
    parser.add_argument('-d', '--distances',
                        help='CSV of pairwise haplotype distances')
    parser.add_argument('-g', '--gfa', help='GFA file with path nodes')
    parser.add_argument('-l', '--log-dir',
                        help='Directory with haplotype sampling logs')
    parser.add_argument('-r', '--reads-dir',
                        help='Directory with TSVs of truth read positions')
    parser.add_argument('-a', '--aln-dir',
                        help='Directory with TSVs for read alignments')
    return parser.parse_args()

def read_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read the cenhap assignment table.
    
    Assuming the first line is a header,
    reads haplotype,cenhap lines into 
    {haplotype : cenhap} dictionary.
    """
    cenhap_table = dict()
    with open(cenhap_file) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            cenhap_table[parts[0]] = parts[1]
    return cenhap_table

def read_distances(distances_file: str) -> Dict[str, Dict[str, float]]:
    """Read the pairwise distance file.

    Three-column CSV with hap1,hap2,dist
    into symmetrical dictionary of
    {hap1 : {hap2 : dist}}.
    """
    dist_matrix = dict()
    with open(distances_file) as file:
        for line in file:
            parts = line.strip().split(',')
            hap1, hap2, dist = parts[0], parts[1], float(parts[2])

            # In case these haplotypes haven't been seen before
            if hap1 not in dist_matrix:
                dist_matrix[hap1] = dict()
            if hap2 not in dist_matrix:
                dist_matrix[hap2] = dict()
            
            dist_matrix[hap1][hap2] = dist
            dist_matrix[hap2][hap1] = dist
    return dist_matrix

def parse_node_list(node_str: str) -> Set[int]:
    """Convert comma-separated list of node IDs to set of normalized IDs."""
    return {int(node.rstrip('+-')) for node 
            in node_str.strip().split(',') if node}

def find_private_nodes(gfa_file: str) -> Dict[str, Set[int]]:
    """Use GFA file to find nodes private to each haplotype.
    
    Reads each path's nodes from the P lines,
    and returns a dictionary with {path : {nodes}}
    for nodes appearing in only one path.
    """

    path_nodes = dict()
    node_counts = dict()

    with open(gfa_file) as f:
        for line in f:
            if not line.startswith('P'):
                continue

            parts = line.rstrip().split()
            # Extract contig name from <sample ID>#<hap num>#<contig>
            path_name = parts[1].split('#')[-1]
            nodes = parse_node_list(parts[2])

            path_nodes[path_name] = nodes
            for node in nodes:
                # In case this node hasn't been seen before
                if node not in node_counts:
                    node_counts[node] = 0
                node_counts[node] += 1

    # Subset each paths' nodes to only the unique ones
    return {path_name : {n for n in nodes if node_counts[n] == 1}
            for path_name, nodes in path_nodes.items()}

def find_guesses(log_file: str) -> Tuple[List[str], str]:
    """Look up the guessed # of sampled haplotypes & cenhap.
    
    Pulls # of sampled haplotypes and guessed cenhap from line
    "Best guess: use <n> haplotypes & sample cenhap = <cenhap>".
    Also pulls specific sampled haplotypes via
    "Selected haplotype <name> with score <score>" lines
    which must appear before the best-guess line.

    Returns a list of the haplotypes used, in order, along with
    the cenhap guessed for the source haplotype.
    """

    sampled_haps = []
    with open(log_file) as file:
        for line in file:
            if line.startswith('Selected haplotype'):
                sampled_haps.append(line.split()[2])
            elif line.startswith('Best guess'):
                parts = line.strip().split()
                n_haps = int(parts[3])
                return (sampled_haps[:n_haps], parts[-1])
            
    return (None, None)

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
            if len(parts) < 2:
                print(read_tsv)
                return None
            truth_nodes[parts[0]] = parse_node_list(parts[1])

    return truth_nodes

def calc_aln_stats(truth_nodes: Dict[str, Set[int]], private_nodes: Set[int],
                   aln_tsv_file: str, is_native: bool) -> Tuple[float, float]:
    """Calculate mean identity & correctness for an alignment set.
    
    Takes in an alignment TSV from vg filter --tsv-out "name;identity;nodes"
    and uses truth/private node lists to get correctness stats.
    Outputs tuple of (mean identity, mean correctness).

    is_native indicates whether these alignments were allowed to use the
    native haplotype; if not, then no penalty is applied when nodes
    private to the haplotype path are missed. See "Correctness calculations".
    """

    identity_scores = []
    correctness_scores = []

    with open(aln_tsv_file) as file:
        file.readline() # header
        for line in file:
            parts = line.strip().split()
            identity_scores.append(float(parts[1]))

            # Calculate correctness of this node list
            if len(parts) == 2 and parts[1] == '0':
                # Unmapped read has 0 correctness
                correctness_scores.append(0)
            elif len(parts) == 3:
                aln_nodes = parse_node_list(parts[2])
                truth = truth_nodes[parts[0]]

                if not is_native:
                    # Don't require alignments to private nodes
                    truth -= private_nodes
                if not truth:
                    # All truth nodes have been thrown out
                    continue
                
                overlap = truth & aln_nodes
                correctness_scores.append(100 * len(overlap) / len(truth))
            else:
                raise ValueError(f'No nodes in line from {aln_tsv_file}')

    if not identity_scores:
        print(aln_tsv_file)
        return None, None
    mean_identity = sum(identity_scores) / len(identity_scores)
    mean_correctness = sum(correctness_scores) / len(correctness_scores)
    return mean_identity, mean_correctness

def write_data(cenhap_table: Dict[str, str], 
               dist_matrix: Dict[str, Dict[str, float]],
               private_nodes: Dict[str, Set[str]]) -> None:
    """Write the output TSV (see file docstring)."""
    # Write header
    column_titles = ['Path name', 'Truth cenhap', 'Guessed cenhap',
                     '# haplotypes sampled', 'Minimum graph distance',
                     'Minimum sampled distance']
    for aln_group in ALN_SUFFIXES.keys():
        column_titles += [f'{aln_group} identity', f'{aln_group} correctness']
    print('\t'.join(column_titles))

    for path_name, truth_cenhap in cenhap_table.items():
        items_to_write = [path_name, truth_cenhap]
        
        guess_file = os.path.join(args.log_dir, f'{path_name}.guess.log')
        if not os.path.exists(guess_file):
            # No haplotype sampling occurred; skip
            continue

        # Add guesses from logfile
        sampled, guess_cenhap = find_guesses(guess_file)
        if not sampled:
            # No haplotype sampling occurred; skip
            continue
        items_to_write += [guess_cenhap, len(sampled)]

        # Look up distances
        dist_row = dist_matrix[path_name]
        # Min distance to any haplotype
        true_closest_hap = min(dist_row, key=dist_row.get)
        # Min distance to a sampled haplotype
        closest_sampled_hap = min(sampled, key=dist_row.get)
        items_to_write += [dist_row[true_closest_hap], 
                           dist_row[closest_sampled_hap]]

        truth_node_file = os.path.join(args.reads_dir, 
                                       f'real_{path_name}.chr12.hifi.tsv')
        if not os.path.exists(truth_node_file):
            continue
        truth_nodes = load_truth_nodes(truth_node_file)
        if truth_nodes is None:
            continue
        for aln_group, tsv_suffix in ALN_SUFFIXES.items():
            # Construct appropriate inputs
            aln_tsv_file = f'{args.aln_dir}/{path_name}.{tsv_suffix}'
            is_native = aln_group.startswith('Native')
            # Add stats to the list of values to write
            identity, correctness = calc_aln_stats(
                truth_nodes, private_nodes[path_name], aln_tsv_file, is_native)
            if identity is None:
                continue
            items_to_write += [identity, correctness]

        print('\t'.join(str(item) for item in items_to_write))

if __name__ == '__main__':
    args = parse_args()
    # Load cross-sample files first
    cenhap_table = read_cenhap_table(args.cenhap_table)
    dist_matrix = read_distances(args.distances)
    private_nodes = find_private_nodes(args.gfa)
    # Write by-sample data
    write_data(cenhap_table, dist_matrix, private_nodes)