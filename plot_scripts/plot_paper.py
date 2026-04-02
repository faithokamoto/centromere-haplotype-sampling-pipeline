#!/usr/bin/env python3

import argparse
import csv
import os
import re

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn
import pandas as pd

CHROM_ORDER = ['chr4', 'chr6', 'chr9', 'chr10', 'chr11', 'chr12', 'chr17']

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-tsv', help='TSV with data to plot')
    parser.add_argument('-o', '--output-dir', help='Directory to write to')
    return parser.parse_args()

def read_tsv(path: str):
    """Read in raw data from TSV."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames

        for r in reader:
            row = dict(r)
            for k in r:
                if k.endswith('cenhap'):
                    row[k] = int(r[k])
                else:
                    try:
                        row[k] = float(r[k])
                    except:
                        pass
            rows.append(row)

    return rows, header

def group_by(rows: dict, col: str) -> dict:
    """Organize rows by their truth cenhap"""
    groups = dict()
    for r in rows:
        cenhap = r[col]
        if not cenhap in groups:
            groups[cenhap] = list()
        groups[cenhap].append(r)
    return groups

def collect_metric(groups: dict, cenhaps: list, columns: list) -> dict:
    """Get raw data for a particular cenhap/col combo."""
    data = {col: [] for col in columns}

    for col in columns:
        for ch in cenhaps:
            vals = [r[col] for r in groups[ch] if isinstance(r[col], float)]
            data[col].append(vals)

    return data

def _chrom_sort_key(val):
    """
    Sorting key for chromosome labels like:
    chr1, chr2, ..., chr22, chrX, chrY, chrM
    """
    if val is None:
        return (99, "")

    s = str(val).lower()

    match = re.match(r"chr(\d+|x|y|m)$", s)
    if match:
        token = match.group(1)
        if token.isdigit():
            return (0, int(token))  # numeric chromosomes first
        special_map = {"x": 23, "y": 24, "m": 25}
        return (0, special_map.get(token, 99))

    # fallback: non-standard labels go last, alphabetical
    return (1, s)

def grouped_violin(ax: plt.Axes, groups: list, data_dict: dict,
                   columns: list, colors: list) -> None:
    """Grouped violin plot with stripplot overlay and chromosome sorting."""

    # ---- sort groups ----
    try:
        sorted_groups = sorted(groups, key=_chrom_sort_key)
    except Exception:
        sorted_groups = groups

    group_index_map = {g: i for i, g in enumerate(groups)}

    width = 0.8
    n = len(columns)
    offsets = np.linspace(-width / 2, width / 2, n)

    for i, col in enumerate(columns):
        reordered_data = []

        for g in sorted_groups:
            original_idx = group_index_map[g]
            values = data_dict[col][original_idx]
            reordered_data.append(values)

        pos = np.arange(len(sorted_groups)) + offsets[i]

        # ---- violin ----
        parts = ax.violinplot(
            reordered_data,
            positions=pos,
            widths=width / n,
            showmeans=True,
            showextrema=False,
        )

        for pc in parts['bodies']:
            pc.set_facecolor(colors[i])
            pc.set_edgecolor('black')
            pc.set_alpha(0.6)

        parts['cmeans'].set_color('black')
        parts['cmeans'].set_linewidth(1.5)

        # ---- stripplot overlay (manual for alignment) ----
        for j, values in enumerate(reordered_data):
            x_center = pos[j]

            # jitter around the violin center
            jitter = (np.random.rand(len(values)) - 0.5) * (width / n * 0.6)
            x_vals = np.full(len(values), x_center) + jitter

            ax.scatter(
                x_vals,
                values,
                color=colors[i],
                alpha=0.6,
                s=10,
                edgecolors='none'
            )

    # ---- axes ----
    ax.set_xticks(range(len(sorted_groups)))
    ax.set_xticklabels(sorted_groups, rotation=45, ha='right')
    ax.set_xlim(-0.5, len(sorted_groups) - 0.5)

def make_heatmap(ax: plt.Axes, rows: dict):
    """Plot heatmap of cenhap typing accuracy."""
    # UNION of all labels → ensures square matrix with same axes
    labels = sorted({r['Truth cenhap'] for r in rows} | {r['Guessed cenhap'] for r in rows})
    index = {v: i for i, v in enumerate(labels)}

    counts = np.zeros((len(labels), len(labels)))

    correct = 0
    total = 0
    for r in rows:
        counts[index[r['Truth cenhap']], index[r['Guessed cenhap']]] += 1
        if r['Truth cenhap'] == r['Guessed cenhap']:
            correct += 1
        total += 1

    print(f'Accuracy: {(100 * correct / total):2f}')

    # Normalize row-wise
    norm = counts.copy()
    for i in range(len(labels)):
        s = counts[i].sum()
        if s > 0:
            norm[i] /= s

    im = ax.imshow(norm, cmap='viridis', vmin=0, vmax=1)

    # annotate counts
    for i in range(len(labels)):
        for j in range(len(labels)):
            if counts[i, j] > 0:
                ax.text(
                    j, i, int(counts[i, j]),
                    ha='center', va='center',
                    backgroundcolor='white', color='black'
                )

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)

    ax.set_xlabel('Guessed cenhap')
    ax.set_ylabel('True cenhap')

    ax.set_box_aspect(1)

    return im

def panel_letter(ax: plt.Axes, letter: str):
    """Add a letter label to a subpanel."""
    ax.text(-0.05, 1.05, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold")

def legend_from_columns(fig: plt.Figure, columns: list, colors: list):
    """Generate a legend for a figure."""
    handles = []
    for c, col in zip(colors, columns):
        # Remove " sim identity" etc. from end of column name
        basename = ' '.join(col.split()[:-2])
        handles.append(mpatches.Patch(color=c, label=basename))
    cols = (len(handles) + 1) // 2
    fig.legend(handles=handles, loc='lower right', ncol=cols, frameon=False)

def filter_columns(header, col_type: str, real: bool,
                   include_all_giraffe: bool = False) -> list:
    """Filter out columns which won't be in the plot."""
    cols = []
    for h in header:
        if not col_type in h:
            continue
        if (real and 'sim' in h) or (not real and 'real' in h):
            continue

        if 'minimap2' in h:
            cols.append(h)
        elif 'giraffe' in h:
            if include_all_giraffe:
                cols.append(h)
            elif 'Sampled giraffe' in h:
                cols.append(h)

    return cols

def fig5(rows, header, outdir):

    filt = [r for r in rows if r["Minimum graph distance"] < 0.2]

    groups = group_by(filt, 'Chromosome')
    chroms = sorted(groups.keys())

    id_cols = filter_columns(header, 'identity', True, False)
    cor_cols = filter_columns(header, 'correctness', True, False)

    id_data = collect_metric(groups, chroms, id_cols)
    cor_data = collect_metric(groups, chroms, cor_cols)

    colors = plt.cm.tab10.colors[:len(id_cols)]

    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(20, 10))

    correct_swarms(axs[0], rows)

    grouped_violin(axs[1], chroms, id_data, id_cols, colors)
    axs[1].set_ylabel("Identity")

    grouped_violin(axs[2], chroms, cor_data, cor_cols, colors)
    axs[2].set_ylabel("Correctness")
    axs[2].set_xlabel("Truth cenhap")

    panel_letter(axs[0], "a")
    panel_letter(axs[1], "b")
    panel_letter(axs[2], "c")

    legend_from_columns(fig, id_cols, colors)

    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(os.path.join(outdir, "fig5.png"), dpi=300)
    plt.close(fig)

def sup_haploid_heatmaps(rows, outdir):
    groups = group_by([r for r in rows], 'Chromosome')

    fig, axs = plt.subplots(2, 4, figsize=(30, 15))

    for i, chrom in enumerate(CHROM_ORDER):
        cur_row = i // 4
        cur_col = i % 4
        make_heatmap(axs[cur_row][cur_col], groups[chrom])
        axs[cur_row][cur_col].set_title(chrom)
        if cur_row == 0:
            axs[cur_row][cur_col].set_ylabel("")

    axs[1][3].set_axis_off()

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "sup_haploids.png"), dpi=300)
    plt.close(fig)


def sup_violin(rows, header, outdir, real, name, dist_filter):
    subset = [r for r in rows if dist_filter(r["Minimum graph distance"])]
    chrom_groups = group_by(subset, 'Chromosome')

    colors = plt.cm.tab20.colors[:8]
    colors = colors[:5] + colors[6:]

    fig, axs = plt.subplots(nrows=4, ncols=len(chrom_groups), figsize=(80, 24))
    
    for i, chrom in enumerate(sorted(chrom_groups.keys())):
        cenhap_groups = group_by(chrom_groups[chrom], 'Truth cenhap')
        cenhaps = sorted(cenhap_groups.keys())

        id_cols = filter_columns(header, 'identity', real, True)
        cor_cols = filter_columns(header, 'correctness', real, True)
        runtime_cols = filter_columns(header, 'runtime', real, True)
        memory_cols = filter_columns(header, 'memory', real, True)

        id_data = collect_metric(cenhap_groups, cenhaps, id_cols)
        cor_data = collect_metric(cenhap_groups, cenhaps, cor_cols)
        runtime_data = collect_metric(cenhap_groups, cenhaps, runtime_cols)
        memory_data = collect_metric(cenhap_groups, cenhaps, memory_cols)

        grouped_violin(axs[0][i], cenhaps, runtime_data, runtime_cols, colors)
        grouped_violin(axs[1][i], cenhaps, memory_data, memory_cols, colors)
        grouped_violin(axs[2][i], cenhaps, id_data, id_cols, colors)
        grouped_violin(axs[3][i], cenhaps, cor_data, cor_cols, colors)

        axs[0][i].set_title(chrom)
        axs[3][i].set_xlabel('Truth cenhap')

    
    axs[0][0].set_ylabel("Runtime", fontsize=13)
    axs[1][0].set_ylabel("Memory", fontsize=13)
    axs[2][0].set_ylabel("Identity", fontsize=13)
    axs[3][0].set_ylabel("Correctness", fontsize=13)
    panel_letter(axs[0][0], "a")
    panel_letter(axs[1][0], "b")
    panel_letter(axs[2][0], "c")
    panel_letter(axs[3][0], "d")

    legend_from_columns(fig, id_cols, colors)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(outdir, name), dpi=300)
    plt.close(fig)

def sup_dists(rows, outdir):
    plt.figure(figsize=(6, 6))

    plt.scatter([r["Minimum graph distance"] for r in rows],
                [r["Minimum sampled distance"] for r in rows],
                s=12, alpha=0.7, color='black')

    plt.axline((0, 0), slope=1)

    plt.xlabel("Minimum graph distance", fontsize=13)
    plt.ylabel("Minimum sampled distance", fontsize=13)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sup_dists.png"), dpi=300)
    plt.close()

def correct_swarms(ax, rows):

    # Convert to DataFrame
    df = pd.DataFrame(rows)

    # Drop rows missing required values
    df = df.dropna(subset=[
        "Chromosome",
        "Minimum graph distance",
        "Truth cenhap",
        "Guessed cenhap"
    ])

    # Add correctness column
    df["Correctness"] = df.apply(
        lambda r: "Correct" if r["Truth cenhap"] == r["Guessed cenhap"] else "Incorrect",
        axis=1
    )

    # Sort chromosomes for nicer plotting (optional)
    try:
        df["Chromosome"] = pd.Categorical(
            df["Chromosome"],
            categories=sorted(df["Chromosome"].unique(), key=_chrom_sort_key),
            ordered=True
        )
    except Exception:
        pass

    seaborn.swarmplot(
        data=df,
        x="Chromosome",
        y="Minimum graph distance",
        hue="Correctness",
        dodge=True,   # <-- this gives two side-by-side swarms per chromosome
        alpha=0.7,
        size=3,
        ax=ax
    )

    ax.set_ylabel("All-pairs distance to neighbor")

def sup_correctness(rows, outdir):
    # Identify all correctness column pairs
    keys = rows[0].keys()
    pairs = []
    for k in keys:
        if "real correctness" in k:
            sim_k = k.replace("real correctness", "sim correctness")
            if sim_k in keys:
                pairs.append((k, sim_k))

    # Color map
    cmap = plt.cm.tab10
    colors = {pair: cmap(i) for i, pair in enumerate(pairs)}

    plt.figure(figsize=(6, 6))

    # Plot each pair
    for (real_k, sim_k), color in colors.items():
        x = [r[sim_k] for r in rows if r[sim_k] is not None and r[real_k] is not None]
        y = [r[real_k] for r in rows if r[sim_k] is not None and r[real_k] is not None]

        label = real_k.replace(" real correctness", "")
        plt.scatter(x, y, s=12, alpha=0.7, color=color, label=label)

    # y = x reference line
    plt.axline((0, 0), slope=1)

    plt.xlabel("Sim correctness", fontsize=13)
    plt.ylabel("Real correctness", fontsize=13)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sup_correctness.png"), dpi=300)
    plt.close()

if __name__ == '__main__':
    args = parse_args()

    rows, header = read_tsv(args.input_tsv)

    fig5(rows, header, args.output_dir)