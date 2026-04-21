#!/usr/bin/env python3
"""Make a detailed alignment comparison supplementary figure.

This figure focuses on the chr10 haplotype for HG01106.1.

====Plot description ====

Rows are vs. CHM13, vs. HG01891.2, vs. a personalized graph, and vs. self
Cols are identity for reads hailing from a location, and depth on linear ref.

==== How I ran this ====

From within /private/home/fokamoto/centromere-haplotype-sampling-pipeline

./hap_focus/alignment_details_fig.py \
    -r /private/groups/patenlab/fokamoto/centrolign/to_align/chr10.HG01106.1.real.truth.sam \
    -c chr10 -H HG01106.1 -N HG01891.2 \
    -a /private/groups/patenlab/fokamoto/centrolign/alignments/haploid \
    -g /private/groups/patenlab/fokamoto/centrolign/graph/haploid \
    -o plot_outputs/aln_details
"""

# ==== file setup ====

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Dict, List, Tuple # Type hinting

import matplotlib.font_manager as fm # Force Arial
import matplotlib.patches as mplpatches # Rectangles
import matplotlib.pyplot as plt # Basic plotting

BIN_SIZE = 10_000
"""How wide are bins for the barplots? 1-->BIN_SIZE is first bin."""

# Use custom TTF for font
arial = fm.FontEntry(fname='./input_data/arial.ttf', name='Arial')
fm.fontManager.ttflist.insert(0, arial)
plt.rcParams['font.family'] = arial.name

# Total figure
figure_size = (13, 3.15) # (width, height)
figure_dpi = 600
# All panels have the same dimensions
panel_width = 3
panel_height = 0.5
# Panel locations
chm13_id_loc = (1.75, 2.3) # (left, bottom)
neighbor_id_loc = (1.75, 1.6)
sampled_id_loc = (1.75, 0.9)
self_id_loc = (1.75, 0.2)

chm13_depth_loc = (5.25, 2.3)
neighbor_depth_loc = (5.25, 1.6)
self_depth_loc = (5.25, 0.2)

hap1_depth_loc = (8.75, 2.3)
hap2_depth_loc = (8.75, 1.6)
hap3_depth_loc = (8.75, 0.9)
hap4_depth_loc = (8.75, 0.2)

# ==== general helpers ====

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--truth-reads-sam', required=True,
                        help='SAM file with original read truth positions')
    parser.add_argument('-c', '--chromosome', required=True,
                        help='Name of chromosome under focus')
    parser.add_argument('-H', '--haplotype-name', required=True,
                        help='Name of haplotype under focus')
    parser.add_argument('-N', '--neighbor-name', required=True,
                        help='Name of neighbor haplotype')
    parser.add_argument('-a', '--aln-dir', required=True,
                        help='Directory with SAMs/TSVs from alignments')
    parser.add_argument('-g', '--graph-dir', required=True,
                        help='Directory with FASTAs/graphs for haplotypes')
    parser.add_argument('-o', '--output-prefix', required=True,
                        help='Prefix for svg/png outputs')
    return parser.parse_args()

def coord_to_bin(coord: int) -> int:
    """Convert raw coordinate to a BIN_SIZE bin start."""
    # Round down, with all of 1-->BIN_SIZE ending up at 0
    bin_index = int((coord - 1) // BIN_SIZE)
    # Convert back to a bin start
    return (bin_index * BIN_SIZE) + 1

# ==== data collection functions ====

def get_hap_len(fasta_file: str) -> int:
    """Calculate length of a haplotype from its FASTA file.
    
    Assumes that only one sequence exists in the FASTA file.
    """
    total_len = 0
    with open(fasta_file) as file:
        file.readline() # Skip '>' name line
        for line in file:
            total_len += len(line.strip())
    return total_len

def read_truth_pos(sam_file: str) -> Dict[str, int]:
    """Read a SAM file into a {read name : start coordinate} dict."""
    read_starts = dict()
    with open(sam_file) as file:
        for line in file:
            parts = line.split('\t')
            read_starts[parts[0]] = int(parts[3])
    return read_starts

def read_coverage(depth_tsv: str, hap_len: int) -> Dict[int, int]:
    """Read a TSV from `samtools depth` into binned average coverage values."""
    bin_cov = dict()
    cur_total_cov = 0
    with open(depth_tsv) as file:
        for line in file:
            parts = line.split('\t')
            pos = int(parts[1])
            cur_total_cov += int(parts[2])
            if pos % BIN_SIZE == 0:
                # We've reached the end of a bin; save average
                bin_cov[coord_to_bin(pos)] = cur_total_cov / BIN_SIZE
                cur_total_cov = 0
        last_bin = coord_to_bin(hap_len)
        last_bin_length = hap_len + 1 - last_bin
        bin_cov[last_bin] = cur_total_cov / last_bin_length
    return bin_cov

def read_identities(aln_tsv: str) -> Dict[str, float]:
    """Read a TSV with name/identity into a {read name : identity} dict."""
    read_id_vals = dict()
    with open(aln_tsv) as file:
        file.readline() # Skip header line
        for line in file:
            parts = line.split('\t')
            read_id_vals[parts[0]] = float(parts[1])
    return read_id_vals

def read_sampled_haps(log_file: str) -> List[str]:
    """Read sampled haplotype names in order.
    
    Pulls specific sampled haplotypes via
        Selected haplotype <name> with score <score>
    """
    sampled_haps = []
    with open(log_file) as file:
        for line in file:
            parts = line.strip().split()
            if line.startswith('Selected haplotype'):
                hap_name = parts[2]
                sampled_haps.append(hap_name)
            elif line.startswith('Best guess'):
                n_haps = int(parts[3])
                return sampled_haps[:n_haps]

def read_walk_nodes(gfa_file: str) -> Dict[int, Tuple[List[int], int]]:
    """Read new nodes along walks into an ordered list.

    Keys are the haplotype number from column 3.
    Values are nodes in order along that walk.
    
    Assumes only '>' characters separate nodes, i.e. nothing backwards.

    Goes through the nodes and adds a negative if it's been seen
    in an earlier-numbered haplotype.

    Returns {hap num : ([nodes], path len)}
    """
    walk_nodes = dict()
    with open(gfa_file) as file:
        for line in file:
            if line.startswith('W'):
                parts = line.split()
                hap_num = int(parts[2])
                hap_len = int(parts[5])
                # Extract node IDs from >1>3>274408>etc.
                nodes = [int(node_id) for node_id in parts[6].split('>')[1:]]
                walk_nodes[hap_num] = (nodes, hap_len)
    
    # Negate previously-seen nodes
    seen_nodes = set()
    for hap_num in sorted(walk_nodes.keys()):
        for i, node in enumerate(walk_nodes[hap_num][0]):
            if node in seen_nodes:
                walk_nodes[hap_num][0][i] = -node
            else:
                seen_nodes.add(node)
    
    return walk_nodes  

def read_pos_depths(depth_tsv: str) -> Dict[int, List[int]]:
    """Read a vg depths --as-table TSV into a {node : {offset : depth}} dict."""
    pos_depths = dict()
    with open(depth_tsv) as file:
        file.readline() # Skip header line
        for line in file:
            parts = line.split('\t')
            node_id = int(parts[1])
            if not node_id in pos_depths:
                pos_depths[node_id] = dict()
            pos_depths[node_id][int(parts[2])] = int(parts[3])
    for node_id, offset_depths in pos_depths.items():
        dep_list = [offset_depths[off] for off in sorted(offset_depths.keys())]
        pos_depths[node_id] = dep_list
    return pos_depths

def bin_walk_depths(pos_depths: Dict[int, List[int]], 
                    hap_nodes: List[int]) -> Dict[int, int]:
    """Get average novel depths within bins along a walk."""
    bin_cov = dict()
    cur_total_cov = 0
    cur_pos = 1
    for node in hap_nodes:
        abs_node = node if node > 0 else -node
        for depth in pos_depths[abs_node]:
            if node > 0:
                cur_total_cov += depth
            if cur_pos % BIN_SIZE == 0:
                # We've reached the end of a bin; save average
                bin_cov[cur_pos] = cur_total_cov / BIN_SIZE
                cur_total_cov = 0
            cur_pos += 1
    last_bin = coord_to_bin(cur_pos - 1)
    last_bin_length = cur_pos - last_bin
    bin_cov[last_bin] = cur_total_cov / last_bin_length
    return bin_cov

# ==== plotting functions ====

def set_up_panel(figsize: Tuple[float, float],
                 loc: Tuple[float, float]) -> plt.Axes:
    """Create a panel at (bottom, left)."""
    panel = plt.axes([loc[0] / figsize[0], loc[1] / figsize[1],
                     panel_width / figsize[0], panel_height / figsize[1]])
    panel.set_xticks([])
    panel.set_yticks([])
    return panel

def row_label(panel: plt.Axes, label: str, left: bool) -> None:
    """Add a row label to the side of a column."""
    x = -0.5 if left else 1.1
    panel.text(x, 0.5, label, transform=panel.transAxes,
               horizontalalignment='left', verticalalignment='center')

def plot_identity_barchart(panel: plt.Axes, read_id_vals: Dict[str, float],
                           read_starts: Dict[str, int], hap_len: int) -> None:
    """Plot binned average identity barchart."""
    # Collect all identity values in each bin
    bins = {bin_start + 1: [] for bin_start 
            in range(0, coord_to_bin(hap_len) + 1, BIN_SIZE)}
    for name, id_val in read_id_vals.items():
        cur_bin = coord_to_bin(read_starts[name])
        bins[cur_bin].append(id_val)

    for bin_start, id_vals in bins.items():
        if id_vals:
            avg_id = sum(id_vals) / len(id_vals)
            panel.add_patch(mplpatches.Rectangle(
                    (bin_start, 0), 
                    BIN_SIZE, avg_id,
                    facecolor='black'
                ))

    # Interesting identity range is between 0.9 and 1
    panel.set_xlim(1, hap_len)
    panel.set_ylim(0.9, 1)
    panel.set_yticks([0.9, 1], ['.9', '1'])

def plot_depth_barchart(panel: plt.Axes, bin_cov: Dict[int, int],
                        hap_len: int) -> None:
    """Plot binned average depth barchart."""
    max_cov = 0
    for bin_start, cur_cov in bin_cov.items():
        max_cov = max(max_cov, cur_cov)
        if cur_cov:
            panel.add_patch(mplpatches.Rectangle(
                    (bin_start, 0), 
                    BIN_SIZE, cur_cov, 
                    facecolor='black'
                ))

    panel.set_xlim(1, hap_len)
    panel.set_ylim(0, max_cov)

    # Floor to one decimal place
    top_tick = int(max_cov * 10) / 10
    panel.set_yticks([0, top_tick])

if __name__ == '__main__':
    args = parse_args()
    hap_prefix = f'{args.chromosome}.{args.haplotype_name}'
    aln_prefix = os.path.join(args.aln_dir, hap_prefix)
    graph_prefix = os.path.join(args.graph_dir, hap_prefix)

    # Get data
    truth_pos = read_truth_pos(args.truth_reads_sam)

    chm13_len = get_hap_len(
        os.path.join(args.graph_dir, f'{args.chromosome}.CHM13.fasta'))
    neighbor_len = get_hap_len(
        os.path.join(args.graph_dir, 
                     f'{args.chromosome}.{args.neighbor_name}.fasta'))
    self_len = get_hap_len(f'{graph_prefix}.fasta')

    chm13_depth = read_coverage(f'{aln_prefix}.CHM13.depth.tsv', chm13_len)
    neighbor_depth = read_coverage(f'{aln_prefix}.neighbor.depth.tsv', 
                                   neighbor_len)
    self_depth = read_coverage(f'{aln_prefix}.own_hap.depth.tsv', self_len)

    chm13_id = read_identities(f'{aln_prefix}.CHM13.real.minimap2.tsv')
    neighbor_id = read_identities(f'{aln_prefix}.neighbor.real.minimap2.tsv')
    sampled_id = read_identities(f'{aln_prefix}.sampled.real.giraffe.tsv')
    self_id = read_identities(f'{aln_prefix}.own_hap.real.minimap2.tsv')
    
    sampled_hap_names = read_sampled_haps(f'{graph_prefix}.guess.real.log')

    walk_nodes = read_walk_nodes(f'{graph_prefix}.sampled.real.gfa')
    pos_depths = read_pos_depths(f'{aln_prefix}.sampled.real.giraffe.pos.depth')
    hap1_bins = bin_walk_depths(pos_depths, walk_nodes[1][0])
    hap2_bins = bin_walk_depths(pos_depths, walk_nodes[2][0])
    hap3_bins = bin_walk_depths(pos_depths, walk_nodes[3][0])
    hap4_bins = bin_walk_depths(pos_depths, walk_nodes[4][0])

    # Make panels
    fig = plt.figure(figsize=figure_size, dpi=figure_dpi)

    chm13_id_panel = set_up_panel(figure_size, chm13_id_loc)
    chm13_depth_panel = set_up_panel(figure_size, chm13_depth_loc)
    
    neighbor_id_panel = set_up_panel(figure_size, neighbor_id_loc)
    neighbor_depth_panel = set_up_panel(figure_size, neighbor_depth_loc)

    sampled_id_panel = set_up_panel(figure_size, sampled_id_loc)

    self_id_panel = set_up_panel(figure_size, self_id_loc)
    self_depth_panel = set_up_panel(figure_size, self_depth_loc)

    hap1_depth_panel = set_up_panel(figure_size, hap1_depth_loc)
    hap2_depth_panel = set_up_panel(figure_size, hap2_depth_loc)
    hap3_depth_panel = set_up_panel(figure_size, hap3_depth_loc)
    hap4_depth_panel = set_up_panel(figure_size, hap4_depth_loc)

    # Do plotting
    plot_identity_barchart(chm13_id_panel, chm13_id, truth_pos, self_len)
    plot_identity_barchart(neighbor_id_panel, neighbor_id, truth_pos, self_len)
    plot_identity_barchart(sampled_id_panel, sampled_id, truth_pos, self_len)
    plot_identity_barchart(self_id_panel, self_id, truth_pos, self_len)

    plot_depth_barchart(chm13_depth_panel, chm13_depth, chm13_len)
    plot_depth_barchart(neighbor_depth_panel, neighbor_depth, neighbor_len)
    plot_depth_barchart(self_depth_panel, self_depth, self_len)

    plot_depth_barchart(hap1_depth_panel, hap1_bins, walk_nodes[1][1])
    plot_depth_barchart(hap2_depth_panel, hap2_bins, walk_nodes[2][1])
    plot_depth_barchart(hap3_depth_panel, hap3_bins, walk_nodes[3][1])
    plot_depth_barchart(hap4_depth_panel, hap4_bins, walk_nodes[4][1])

    # Labels for cols 1-2
    row_label(chm13_id_panel, 'CHM13', True)
    row_label(neighbor_id_panel, args.neighbor_name, True)
    row_label(sampled_id_panel, 'Personalized graph', True)
    row_label(self_id_panel, args.haplotype_name, True)

    # Labels for col 3
    row_label(hap1_depth_panel, sampled_hap_names[0], False)
    row_label(hap2_depth_panel, sampled_hap_names[1], False)
    row_label(hap3_depth_panel, sampled_hap_names[2], False)
    row_label(hap4_depth_panel, sampled_hap_names[3], False)

    # Column labels
    chm13_id_panel.set_title('Alignment identity')
    chm13_depth_panel.set_title('Read depth')
    hap1_depth_panel.set_title('Read depth on novel sequence')

    self_id_panel.set_xlabel('Truth position of read')
    self_depth_panel.set_xlabel('Alignment position along linear ref')
    hap4_depth_panel.set_xlabel('Alignment position along sampled hap')

    fig.savefig(f'{args.output_prefix}.png')
    fig.savefig(f'{args.output_prefix}.svg')