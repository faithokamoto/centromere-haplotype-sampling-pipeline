#!/usr/bin/env python3
"""Calculate identity & correctness stats for reads.

== Inputs ==

- GFA file with haplotype paths (--gfa).
    Expects all haplotypes as W lines; ignores all other lines.
    https://gfa-spec.github.io/GFA-spec/GFA1.html
- Truth read positions (<--reads-dir>/<chrom>.<hap>.<realness>.tsv).
    Expects columns for read name and then comma-separated nodes.
    The node list is allowed to end with a comma.
    Also expects a header line but it will be skipped.
- Alignment metadata (<--aln-dir>/<chrom>.<hap>.<ref>.<realness>.<tool>.tsv).
    Expects columns for read name, identity, and comma-separated nodes.
    The node list is allowed to end with a comma.
    Also expects a header line but it will be skipped.
- Alignment logs (<--aln-dir>/<chrom>.<hap>.<ref>.<realness>.<tool>.log).
    Uses /usr/bin/time output to get runtime and memory usage stats.

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
import os # Filesystem interactions
import re # Parsing walks (annoying format)
from typing import Dict, List, Set, Tuple # Type hints

REFS = [('CHM13', 'chm13'), 
        ('Neighbor', 'neighbor'), 
        ('Sampled', 'sampled'), 
        ('Native hap', 'own_hap')]
"""References used in alignments."""
REALNESS = ['real', 'sim']
"""Possible levels of realness."""
ALIGNERS = ['minimap2', 'giraffe']
"""Names of aligners used."""

def parse_args() -> argparse.Namespace:
    """Handle command-line argument parsing."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chrom', help='Chromosome (in filenames)')
    parser.add_argument('-g', '--gfa', help='GFA file with path nodes')
    parser.add_argument('-n', '--hap-name', help='Haplotype reads are from')
    parser.add_argument('-r', '--reads-dir', help='Directory with truth TSVs')
    parser.add_argument('-a', '--aln-dir', help='Directory with aligments')
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

def find_private_nodes(all_nodes: Dict[str, Set[int]], hap_name: str,
                       haps_to_consider: Set[str] = None) -> Set[int]:
    """Find which nodes are private to a particular haplotype.
    
    If haps_to_consider is provided, only consider those haplotypes.
    """
    seen_elsewhere = set()
    if not haps_to_consider:
        # Default to all other haplotypes
        haps_to_consider = set(all_nodes.keys())
    
    haps_to_consider.remove(hap_name)

    for other in all_nodes.keys():
        if other in haps_to_consider:
            seen_elsewhere |= all_nodes[other]
    
    return (all_nodes[hap_name] - seen_elsewhere)

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

def get_guesses(log_file: str) -> List[str]:
    """Look up the sampled haplotypes.
    
    Pulls specific sampled haplotypes via
        Selected haplotype <name> with score <score>
    And grabs the exact number used via
        Best guess: use <n> haplotypes & sample cenhap = <cenhap>

    Returns a list of the haplotypes selected, in order
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

def calc_identity_correctness(truth_nodes: Dict[str, Set[int]], 
                              private_nodes: Set[int],
                              aln_tsv_file: str, 
                              req_private: bool) -> Tuple[float, float]:
    """Calculate mean identity & correctness for an alignment set.
    
    Takes in an alignment TSV from vg filter --tsv-out "name;identity;nodes"
    and uses truth/private node lists to get correctness stats.
    Returns (identity, correctness)

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
                    truth = set(truth_nodes[parts[0]])

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
    mean_identity = sum(identity) / len(identity)
    if correctness:
        mean_correctness = sum(correctness) / len(correctness)
    else:
        mean_correctness = None
    return mean_identity, mean_correctness

def get_runtime_memory(aln_log_file: str) -> Tuple[float, float]:
    """Extract runtime & memory usage stats from logfile.
    
    Expects output from /usr/bin/time in particular using these lines:

    	User time (seconds): X
	    System time (seconds): Y
        ...
	    Maximum resident set size (kbytes): Z

    X + Y is the runtime in seconds, and Z / 10^6 is the memory in GB

    Returns (runtime, memory)
    """

    runtime = 0
    memory = 0

    with open(aln_log_file) as file:
        for line in file:
            if 'User time' in line or 'System time' in line:
                runtime += float(line.strip().split()[-1])
            elif 'Maximum resident set size' in line:
                # Convert from KB to GB
                memory = float(line.strip().split()[-1]) / 1_000_000

    return runtime, memory


def print_aln_stats(truth_nodes: Dict[str, Set[int]], private_nodes: Set[int],
                    aln_prefix: str, req_private: bool) -> None:
    """Print statistics for each alignment run.

    Prints prefix without directory, e.g. chr12.HG00099.1.sampled.real.giraffe
    followed by its stats identity, correctness, runtime, memory
    """

    identity, correctness = calc_identity_correctness(
        truth_nodes, private_nodes, f'{aln_prefix}.tsv', req_private)
    runtime, memory = get_runtime_memory(f'{aln_prefix}.log')
    no_dir_prefix = aln_prefix.split('/')[-1]

    print(f'{no_dir_prefix}: identity {identity} / correctness {correctness} /'
          f' runtime {runtime} / memory {memory}')

if __name__ == '__main__':
    args = parse_args()
    # Needed across alignment files
    all_nodes = read_nodes(args.gfa)
    priv_nodes = find_private_nodes(all_nodes, args.hap_name)
    
    for realness in ['real', 'sim']:
        # Needed across alignment files of a realness
        truth_nodes = load_truth_nodes(os.path.join(args.reads_dir, 
                                f'{args.chrom}.{args.hap_name}.{realness}.tsv'))
        
        for tool in ['minimap2', 'giraffe']:
            for ref in ['CHM13', 'neighbor', 'sampled', 'own_hap']:
                if tool == 'minimap2' and ref == 'sampled':
                    # Minimap2 can't do graph alignments
                    continue

                # Only alignments to a matched assembly must use private nodes
                req_priv = (ref == 'own_hap')
                aln_prefix = os.path.join(args.aln_dir,
                        f'{args.chrom}.{args.hap_name}.{ref}.{realness}.{tool}')
                print_aln_stats(truth_nodes, priv_nodes, aln_prefix, req_priv)