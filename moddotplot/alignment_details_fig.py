#!/usr/bin/env python3
"""Make a detailed alignment comparison supplementary figure.

This figure focuses on the chr10 haplotype for HG01106.1.

====Plot description ====

Panels A (vs. CHM13), B (vs. HG01891.2), and D (vs. itself) are the same
kind of plot.
- The main subpanel has a similarity comparison as output by ModDotPlot,
overlaid with CIGAR string from the Centrolign alignment.
- The top subpanel has a bar plot with average depth in 10kb bins.
- The right-hand subpanel has a bar plot with average identity in 10kb bins.

Panel C still has the same right-hand subpanel, but no top subpanel,
and its main panel is an image of the personalized graph colored by depth.

==== How I ran this ====

From within /private/home/fokamoto/centromere-haplotype-sampling-pipeline

./moddotplot/alignment_details_fig.py \
    -r /private/groups/patenlab/fokamoto/centrolign/to_align/chr10.HG01106.1.real.truth.sam \
    -c chr10 -H HG01106.1 -N HG01891.2 \
    -a /private/groups/patenlab/fokamoto/centrolign/alignments/haploid \
    -f /private/groups/patenlab/fokamoto/centrolign/graph/haploid \
    -m plot_outputs \
    -cc input_data/pairwise_cigar_CHM13.0_HG01106.1.txt \
    -nc input_data/pairwise_cigar_HG01106.1_HG01891.2.txt \
    -o plot_outputs/aln_details
"""

# ==== file setup ====

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Dict, Tuple # Type hinting

import matplotlib.font_manager as fm # Force Arial
import matplotlib.patches as mplpatches # Rectangles
import matplotlib.pyplot as plt # Basic plotting

BIN_SIZE = 10_000
"""How wide are bins for the barplots?"""

# Use custom TTF for font
arial = fm.FontEntry(fname='./input_data/arial.ttf', name='Arial')
fm.fontManager.ttflist.insert(0, arial)
plt.rcParams['font.family'] = arial.name

# Total figure
figure_size = (8.75, 8.75) # (width, height)
figure_dpi = 600
# Panel locations
panelA_top_loc = (0.45, 7.8, 3, 0.5) # (left, bottom, width, height)
panelA_main_loc = (0.45, 4.75, 3, 3)
panelA_right_loc = (3.5, 4.75, 0.5, 3)

panelB_top_loc = (4.75, 7.8, 3, 0.5)
panelB_main_loc = (4.75, 4.75, 3, 3)
panelB_right_loc = (7.8, 4.75, 0.5, 3)

panelC_main_loc = (0.45, 0.45, 3, 3)
panelC_right_loc = (3.5, 0.45, 0.5, 3)

panelD_top_loc = (4.75, 3.5, 3, 0.5)
panelD_main_loc = (4.75, 0.45, 3, 3)
panelD_right_loc = (7.8, 0.45, 0.5, 3)

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
    parser.add_argument('-f', '--fasta-dir', required=True,
                        help='Directory with FASTAs for haplotypes')
    parser.add_argument('-m', '--moddotplot-dir', required=True,
                        help='Directory with ModDotPlot BEDs')
    parser.add_argument('-cc', '--chm13-cigar', required=True,
                        help='Centrolign CIGAR for assembly vs. CHM13')
    parser.add_argument('-nc', '--neighbor-cigar', required=True,
                        help='Centrolign CIGAR for assembly vs. neighbor')
    parser.add_argument('-o', '--output-prefix', required=True,
                        help='Prefix for svg/png outputs')
    return parser.parse_args()

# === data collection functions ====

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

def read_identities(aln_tsv: str) -> Dict[str, float]:
    """Read a TSV with name/identity into a {read name : identity} dict."""
    read_id_vals = dict()
    with open(aln_tsv) as file:
        file.readline() # Skip header line
        for line in file:
            parts = line.split('\t')
            read_id_vals[parts[0]] = float(parts[1])
    return read_id_vals

# ==== plotting functions ====

def set_up_panel(figsize: Tuple[float, float],
                 loc: Tuple[float, float, float, float]) -> plt.Axes:
    """Create a panel at (bottom, left, width, height)."""
    panel = plt.axes([loc[0] / figsize[0], loc[1] / figsize[1],
                     loc[2] / figsize[0], loc[3] / figsize[1]])
    return panel

def panel_letter(panel: plt.Axes, letter: str) -> None:
    """Add a panel letter to the top-left corner."""
    panel.text(-0.15, 1.05, letter, transform=panel.transAxes)

def coord_to_bin(coord: float) -> int:
    """Convert raw coordinate to a BIN_SIZE bin start."""
    return int(coord // BIN_SIZE) * BIN_SIZE

def plot_identity_barchart(panel: plt.Axes, read_id_vals: Dict[str, float],
                           read_starts: Dict[str, int], hap_len: int) -> None:
    """Plot binned average identity sideways barchart."""
    # Interesting identity range is between 0.9 and 1
    panel.set_xlim(0.9, 1)
    panel.set_ylim(0, hap_len)
    panel.set_yticks([])
    # Collect all identity values in each bin
    bins = {bin_start: [] for bin_start in range(0, hap_len, BIN_SIZE)}
    for name, id_val in read_id_vals.items():
        cur_bin = coord_to_bin(read_starts[name])
        bins[cur_bin].append(id_val)

    for bin_start, id_vals in bins.items():
        if id_vals:
            avg_id = sum(id_vals) / len(id_vals)
            panel.add_patch(mplpatches.Rectangle(
                    (0, bin_start), 
                    avg_id, BIN_SIZE, 
                    facecolor='black'
                ))

if __name__ == '__main__':
    args = parse_args()
    hap_prefix = f'{args.chromosome}.{args.haplotype_name}'

    # Get data
    hap_len = get_hap_len(os.path.join(args.fasta_dir, f'{hap_prefix}.fasta'))
    truth_pos = read_truth_pos(args.truth_reads_sam)

    chm13_id_vals = read_identities(
        os.path.join(args.aln_dir, f'{hap_prefix}.CHM13.real.minimap2.tsv'))
    neighbor_id_vals = read_identities(
        os.path.join(args.aln_dir, f'{hap_prefix}.neighbor.real.minimap2.tsv'))
    sampled_id_vals = read_identities(
        os.path.join(args.aln_dir, f'{hap_prefix}.sampled.real.giraffe.tsv'))
    self_id_vals = read_identities(
        os.path.join(args.aln_dir, f'{hap_prefix}.own_hap.real.minimap2.tsv'))
    
    # Make panels
    fig = plt.figure(figsize=figure_size, dpi=figure_dpi)

    panelA_top = set_up_panel(figure_size, panelA_top_loc)
    panelA_main = set_up_panel(figure_size, panelA_main_loc)
    panelA_right = set_up_panel(figure_size, panelA_right_loc)
    
    panelB_top = set_up_panel(figure_size, panelB_top_loc)
    panelB_main = set_up_panel(figure_size, panelB_main_loc)
    panelB_right = set_up_panel(figure_size, panelB_right_loc)
    
    panelC_main = set_up_panel(figure_size, panelC_main_loc)
    panelC_right = set_up_panel(figure_size, panelC_right_loc)
    
    panelD_top = set_up_panel(figure_size, panelD_top_loc)
    panelD_main = set_up_panel(figure_size, panelD_main_loc)
    panelD_right = set_up_panel(figure_size, panelD_right_loc)

    # Do plotting

    plot_identity_barchart(panelA_right, chm13_id_vals, truth_pos, hap_len)
    plot_identity_barchart(panelB_right, neighbor_id_vals, truth_pos, hap_len)
    plot_identity_barchart(panelC_right, sampled_id_vals, truth_pos, hap_len)
    plot_identity_barchart(panelD_right, self_id_vals, truth_pos, hap_len)

    panel_letter(panelA_main, 'a')
    panel_letter(panelB_main, 'b')
    panel_letter(panelC_main, 'c')
    panel_letter(panelD_main, 'd')

    fig.savefig(f'{args.output_prefix}.png')
    fig.savefig(f'{args.output_prefix}.svg')