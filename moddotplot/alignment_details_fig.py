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
    -r /private/groups/patenlab/fokamoto/centrolign/to_align/HG01106.1.real.truth.sam \
    -H HG01106.1 -N HG01891.2 \
    -a /private/groups/patenlab/fokamoto/centrolign/alignments/haploid \
    -m plot_outputs \
    -cc input_data/pairwise_cigar_CHM13.0_HG01106.1.txt \
    -nc input_data/pairwise_cigar_HG01106.1_HG01891.2.txt \
    -o plot_outputs/aln_details
"""

# ==== file setup ====

import argparse # Command-line argument parsing
from typing import Tuple # Type hinting

import matplotlib.font_manager as fm # Force Arial
import matplotlib.patches as mplpatches # Rectangles
import matplotlib.pyplot as plt # Basic plotting

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
    parser.add_argument('-H', '--haplotype-name', required=True,
                        help='Name of haplotype under focus')
    parser.add_argument('-N', '--neighbor-name', required=True,
                        help='Name of neighbor haplotype')
    parser.add_argument('-a', '--aln-dir', required=True,
                        help='Directory with SAMs/TSVs from alignments')
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

if __name__ == '__main__':
    args = parse_args()
    
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

    panel_letter(panelA_main, 'a')
    panel_letter(panelB_main, 'b')
    panel_letter(panelC_main, 'c')
    panel_letter(panelD_main, 'd')

    fig.savefig(f'{args.output_prefix}.png')
    fig.savefig(f'{args.output_prefix}.svg')