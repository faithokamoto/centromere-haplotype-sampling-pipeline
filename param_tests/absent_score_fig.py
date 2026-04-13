#!/usr/bin/env python3
"""Make --absent-score supplementary figure.

==== Panel A ====

Look up all the top-haplotype depth statistics in --depth-log-dir,
then use those to make a grouped histogram calling out how
HG00738.1 is problematic in chr4.

==== Panel B ====

Take a --kmer-count-file with graph-unique k-mer counts,
and plot the distribution among haplotypes.

Assume that the --kmer-count-file lines appear in
increasing count order.

==== Panel C ====

Go through haplotype sampling log files in --sampling-log-dir,
and plot # of times HG00738.1 tricked the typer.

==== How I ran this ====

From within /private/home/fokamoto/centromere-haplotype-sampling-pipeline

./param_tests/absent_score_fig.py \
    -c input_data/chr4.cenhap_predictions.tsv \
    -d /private/groups/patenlab/fokamoto/centrolign/graph/default \
    -k /private/groups/cgl/jlsiren/centrolign/chr4.txt \
    -s log/typing_tests -n HG00738.1 -o plot_outputs/absent_score
"""

# ==== file setup ====

import argparse # Command-line argument parsing
import os # Filesystem interactions
from typing import Dict, List, Tuple # Type hinting

import matplotlib.font_manager as fm # Force Arial
import matplotlib.patches as mplpatches # Rectangles
import matplotlib.pyplot as plt # Basic plotting

# Use custom TTF for font
arial = fm.FontEntry(fname='./param_tests/arial.ttf', name='Arial')
fm.fontManager.ttflist.insert(0, arial)
plt.rcParams['font.family'] = arial.name

DEPTH_SUFFIX = '.guess.log'
"""Suffix to expect for the depth stats logs."""
OTHER_COLOR = '#009E73'
"""Color to use for non-HG00738.1 haplotypes."""
SPECIAL_COLOR = '#E6A024'
"""Color to use for HG00738.1."""
DEPTH_MAX = 30
"""Cap for values in depth histogram."""

# Total figure
figure_size = (10, 2.75) # (width, height)
figure_dpi = 600
# Panel locations
panelA_loc = (0.5, 0.5, 2.5, 2) # (left, bottom, width, height)
panelB_loc = (3.75, 0.5, 2.5, 2)
panelC_loc = (7, 0.5, 2.5, 2)
# Panel limits
panelA_lim = {'x': (-0.5, 30.5), 'y': (0, 16)}
panelB_lim = {'x': (0, 132), 'y': (0, 3700)}
panelC_lim = {'x': (-0.01, 0.81), 'y': (-0.5, 17.5)}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cenhap-table', required=True,
                        help='TSV of truth cenhap assignments')
    parser.add_argument('-d', '--depth-log-dir', required=True,
                        help='Directory with depth stat log files')
    parser.add_argument('-k', '--kmer-count-file', required=True,
                        help='K-mer stats file')
    parser.add_argument('-s', '--sampling-log-dir', required=True,
                        help='Directory with depth stat log files')
    parser.add_argument('-n', '--haplotype-name', required=True,
                        help='Haplotype to call out')
    parser.add_argument('-o', '--output-prefix', required=True,
                        help='Prefix for svg/png outputs')
    return parser.parse_args()

# === data collection functions ====

def read_depth_stats(log_file: str) -> Tuple[str, float]:
    """Read the top haplotype & its private depth from log file.
    
    Look for a line like
        <hap name> has avg depth <depth>
    and return (name, depth)
    """
    with open(log_file) as file:
        for line in file:
            if 'has avg depth' in line:
                parts = line.strip().split()
                top_name = parts[0]
                top_depth = float(parts[4])
                return top_name, top_depth
            
def collect_hap_depths(depth_log_dir: str, depth_log_suffix: str,
                       hap_name: str) -> Tuple[List[int], List[int]]:
    """Catalog depths on top haplotype by special vs. non.

    For all files within the given directory, read the top hap name & depth.
    Return two lists, one with depths for when the top hap == hap_name
    and one with depths on all other top haplotypes.
    """

    special_hap_depths = []
    other_hap_depths = []

    for item in os.listdir(depth_log_dir):
        if item.endswith(depth_log_suffix):
            hap, depth = read_depth_stats(os.path.join(depth_log_dir, item))
            if hap == hap_name:
                special_hap_depths.append(depth)
            else:
                other_hap_depths.append(depth)
    
    return special_hap_depths, other_hap_depths

def read_kmer_counts(log_file: str) -> Tuple[List[int], str]:
    """Read graph-unique k-mer counts on each haplotype.
    
    Look for lines like
        <sample ID>#<hap num>#<hap name>#0: <#> / <total> kmers present (density <density>)
    Save all counts (the <#>) in a list. Also return the name of
    the first haplotype seen.
    """

    kmer_counts = []
    first_hap_name = None
    with open(log_file) as file:
        for line in file:
            if '#0:' in line:
                parts = line.strip().split()
                if first_hap_name is None:
                    first_hap_name = parts[0].split('#')[2]
                kmer_counts.append(int(parts[1]))
    return kmer_counts, first_hap_name

def read_cenhap_table(cenhap_file: str) -> Dict[str, str]:
    """Read a cenhap assignment table as {hap name : truth cenhap}"""
    cenhap_table = dict()
    with open(cenhap_file) as file:
        file.readline()  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            cenhap_table[parts[0]] = parts[1]
    return cenhap_table

def parse_param_test_filename(filename: str) -> Tuple[str, float]:
    """Extract haplotype and --absent-score tested from filename.

    Parse filename of form chr4.<haplotype>.<##>.log
    to return (haplotype, --absent-score)
    """
    parts = filename.split('.')
    # parts = ["chr4", sample_id, hap_num, "##", ".log"]
    score_part = parts[3]
    haplotype = '.'.join(parts[1:3])

    return haplotype, (int(score_part) / 100.0)

def read_first_selected_hap(filepath: str) -> str:
    """Extract first selected haplotype from log file.

    Use lines like
        Selected haplotype <other hap> with score <score>
    """
    with open(filepath) as file:
        for line in file:
            if 'Selected' in line:
                return line.strip().split()[2]
            
def collect_param_test_results(cenhap_table: Dict[str, str], hap_name: str,
                               sampling_log_dir: str) -> Dict[float, int]:
    """Collect correctness results from --absent-score parameter test.
    
    For all files within the given directory,
    parse out their position in the parameter test,
    and use that to figure out if they guessed wrongly
    because of the specific problematic haplotype.

    Return a {--absent-score : # problematic incorrect} dict.
    """
    results = dict()

    for filename in os.listdir(sampling_log_dir):
        testing_hap, absent_score = parse_param_test_filename(filename)
        if not absent_score in results:
            results[absent_score] = 0

        full_filename = os.path.join(sampling_log_dir, filename)
        selected_hap = read_first_selected_hap(full_filename)
        true_cenhap = cenhap_table[testing_hap]
        selected_cenhap = cenhap_table[selected_hap]

        typing_success = (true_cenhap == selected_cenhap)

        if (selected_hap == hap_name) and not typing_success:
            results[absent_score] += 1
    return results

# ==== plotting functions ====

def set_up_panel(figsize: Tuple[float, float],
                 loc: Tuple[float, float, float, float],
                 lim: Dict[str, Tuple[float, float]]) -> plt.Axes:
    """Create a panel at (bottom, left, width, height)."""
    panel = plt.axes([loc[0] / figsize[0], loc[1] / figsize[1],
                     loc[2] / figsize[0], loc[3] / figsize[1]])
    panel.set_xlim(lim['x'])
    panel.set_ylim(lim['y'])
    panel.spines['top'].set_visible(False)
    panel.spines['right'].set_visible(False)
    return panel

def panel_letter(panel: plt.Axes, letter: str) -> None:
    """Add a panel letter to the top-left corner."""
    panel.text(-0.15, 1.05, letter, transform=panel.transAxes)

def depth_to_bin(depth: float) -> float:
    """Convert raw depth to a histogram bin value.
    
    Cap to DPETH_MAX (30) and then round down to an 0.5 increment.
    """
    return int(min(depth, DEPTH_MAX) * 2) / 2

def plot_depth_hist_panel(panel: plt.Axes, special_hap_depths: List[float],
                          other_hap_depths: List[float], hap_name: str) -> None:
    """Plot stacked histogram of depth on top haplotype.
    
    Bottom is non-HG00738.1 haplotypes, and top is HG00738.1.
    Binned in 0.5 increments with >30 values collapsed into 30.
    """

    # Collect bin data
    bins = [val / 2 for val in range(0, DEPTH_MAX * 2 + 1)]
    special_counts = {bin_min : 0 for bin_min in bins}
    other_counts = {bin_min : 0 for bin_min in bins}
    for depth in special_hap_depths:
        special_counts[depth_to_bin(depth)] += 1
    for depth in other_hap_depths:
        other_counts[depth_to_bin(depth)] += 1

    for bin_min in bins:
        if other_counts[bin_min]:
            # Draw bottom of stacked histogram
            panel.add_patch(mplpatches.Rectangle(
                (bin_min - 0.25, 0), 
                0.5, other_counts[bin_min], 
                edgecolor='black', facecolor=OTHER_COLOR
            ))
        if special_counts[bin_min]:
            # Draw top of stacked histogram
            panel.add_patch(mplpatches.Rectangle(
                (bin_min - 0.25, other_counts[bin_min]), 
                0.5, special_counts[bin_min], 
                edgecolor='black', facecolor=SPECIAL_COLOR
            ))
    
    # Manual legend
    panel.add_patch(mplpatches.Rectangle(
        (15, 14), 4, 0.5, edgecolor='black', facecolor=SPECIAL_COLOR
    ))
    panel.add_patch(mplpatches.Rectangle(
        (15, 12), 4, 0.5, edgecolor='black', facecolor=OTHER_COLOR
    ))
    panel.text(20, 14, hap_name)
    panel.text(20, 12, 'Other')

    panel.set_yticks([4, 8, 12, 16])
    panel.set_xlabel('Depth on top haplotype')

def plot_kmer_panel(panel: plt.Axes, kmer_counts: List[int],
                    hap_name: str) -> None:
    """Plot counts on Y with haplotypes sorted along X.
    
    Also call out the minimum k-mer haplotype with a
    larger, specially-colored point.
    """
    panel.plot(range(5, len(kmer_counts) + 5), kmer_counts, color=OTHER_COLOR,
               marker='.', linestyle='None')
    panel.plot(5, kmer_counts[0], color=SPECIAL_COLOR, marker='.')
    panel.text(10, kmer_counts[0] / 2, hap_name)
    panel.set_xticks([])
    panel.set_yticks([1000, 2000, 3000], labels=['1k', '2k', '3k'])
    panel.set_ylabel('k-mers')

def plot_param_search_panel(panel: plt.Axes, error_counts: Dict[float, int],
                            hap_name: str) -> None:
    """Plot errors on Y with --absent-score on X.
    
    Also draw and label a line at x=0.05, where the dots hit y=0.
    """
    absent_scores = sorted(error_counts.keys())
    special_errors = [error_counts[abs_score] for abs_score in absent_scores]
    panel.axvline(x=0.05, color='black', linestyle='--')
    panel.plot(absent_scores, special_errors, color=SPECIAL_COLOR, 
               marker='.', linestyle='None')
    panel.text(0.07, 15, '0.05')

    panel.set_xlabel('--absent-score')
    panel.set_ylabel(f'{hap_name} tricks')

if __name__ == '__main__':
    args = parse_args()

    # Grab all data
    special_hap_depths, other_hap_depths = collect_hap_depths(
        args.depth_log_dir, DEPTH_SUFFIX, args.haplotype_name)
    kmer_counts, min_hap = read_kmer_counts(args.kmer_count_file)
    cenhap_table = read_cenhap_table(args.cenhap_table)
    error_counts = collect_param_test_results(
        cenhap_table, args.haplotype_name, args.sampling_log_dir)
    
    if min_hap != args.haplotype_name:
        raise ValueError(f'--haplotype-name is not the one with least k-mers')
    
    # Make panels
    fig = plt.figure(figsize=figure_size, dpi=figure_dpi)
    panelA = set_up_panel(figure_size, panelA_loc, panelA_lim)
    panelB = set_up_panel(figure_size, panelB_loc, panelB_lim)
    panelC = set_up_panel(figure_size, panelC_loc, panelC_lim)

    plot_depth_hist_panel(panelA, special_hap_depths, other_hap_depths, min_hap)
    plot_kmer_panel(panelB, kmer_counts, min_hap)
    plot_param_search_panel(panelC, error_counts, min_hap)

    panel_letter(panelA, 'a')
    panel_letter(panelB, 'b')
    panel_letter(panelC, 'c')

    fig.savefig(f'{args.output_prefix}.png')
    fig.savefig(f'{args.output_prefix}.svg')