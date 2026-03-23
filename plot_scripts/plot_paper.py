#!/usr/bin/env python3

import argparse
import csv
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

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
                # Try to number-fy the non-cenhap cols
                # (cenhaps need to be raw)
                if not k.endswith('cenhap'):
                    try:
                        row[k] = float(r[k])
                    except:
                        pass
            rows.append(row)

    return rows, header

def group_by_cenhap(rows: dict) -> dict:
    """Organize rows by their truth cenhap"""
    groups = dict()
    for r in rows:
        cenhap = r['Truth cenhap']
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

def grouped_violin(ax: plt.Axes, cenhaps: list, data_dict: dict,
                   columns: list, colors: list) -> None:
    """Make a violin plot grouped by cenhap."""
    width = 0.8
    n = len(columns)
    offsets = np.linspace(-width / 2, width / 2, n)

    for i, col in enumerate(columns):
        pos = np.arange(len(cenhaps)) + offsets[i]

        parts = ax.violinplot(
            data_dict[col],
            positions=pos,
            widths=width / n,
            showmeans=True,
            showextrema=False,
        )

        for pc in parts['bodies']:
            pc.set_facecolor(colors[i])
            pc.set_edgecolor('black')
            pc.set_alpha(0.85)
        
        parts['cmeans'].set_color('black')
        parts['cmeans'].set_linewidth(1.5)

    ax.set_xticks(range(len(cenhaps)))
    ax.set_xticklabels(cenhaps, rotation=45, ha='right')
    ax.set_xlim(-0.5, len(cenhaps) - 0.5)

def make_heatmap(ax: plt.Axes, rows: dict):
    """Plot heatmap of cenhap typing accuracy."""
    truth = sorted({r['Truth cenhap'] for r in rows})
    guess = sorted({r['Guessed cenhap'] for r in rows})

    t_index = {v: i for i, v in enumerate(truth)}
    g_index = {v: i for i, v in enumerate(guess)}

    counts = np.zeros((len(truth), len(guess)))

    for r in rows:
        counts[t_index[r['Truth cenhap']], g_index[r['Guessed cenhap']]] += 1

    norm = counts.copy()
    for i in range(len(truth)):
        s = counts[i].sum()
        if s > 0:
            norm[i] /= s

    im = ax.imshow(norm, cmap='viridis')

    for i in range(len(truth)):
        for j in range(len(guess)):
            if counts[i, j] > 0:
                ax.text(j, i, int(counts[i, j]), ha='center', va='center', 
                        backgroundcolor='white', color='black')

    ax.set_xticks(range(len(guess)))
    ax.set_xticklabels(guess, rotation=45, ha='right')
    ax.set_yticks(range(len(truth)))
    ax.set_yticklabels(truth)

    ax.set_xlabel('Guessed cenhap')
    ax.set_ylabel('Truth cenhap')

    ax.set_box_aspect(1)  # make square

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
    fig.legend(handles=handles, loc='upper right', ncol=cols, frameon=False)

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

    groups = group_by_cenhap(filt)
    cenhaps = sorted(groups.keys())

    id_cols = filter_columns(header, 'identity', True, False)
    cor_cols = filter_columns(header, 'correctness', True, False)

    id_data = collect_metric(groups, cenhaps, id_cols)
    cor_data = collect_metric(groups, cenhaps, cor_cols)

    colors = plt.cm.tab10.colors[:len(id_cols)]

    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.5])

    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 1])

    make_heatmap(axA, rows)

    grouped_violin(axB, cenhaps, id_data, id_cols, colors)
    axB.set_ylabel("Identity")

    grouped_violin(axC, cenhaps, cor_data, cor_cols, colors)
    axC.set_ylabel("Correctness")
    axC.set_xlabel("Truth cenhap")

    panel_letter(axA, "a")
    panel_letter(axB, "b")
    panel_letter(axC, "c")

    legend_from_columns(fig, id_cols, colors)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(outdir, "fig5.png"), dpi=300)
    plt.close(fig)


def sup_violin(rows, header, outdir, real, name, dist_filter):

    subset = [r for r in rows if dist_filter(r["Minimum graph distance"])]

    groups = group_by_cenhap(subset)
    cenhaps = sorted(groups.keys())

    id_cols = filter_columns(header, 'identity', real, True)
    cor_cols = filter_columns(header, 'correctness', real, True)
    runtime_cols = filter_columns(header, 'runtime', real, True)
    memory_cols = filter_columns(header, 'memory', real, True)

    id_data = collect_metric(groups, cenhaps, id_cols)
    cor_data = collect_metric(groups, cenhaps, cor_cols)
    runtime_data = collect_metric(groups, cenhaps, runtime_cols)
    memory_data = collect_metric(groups, cenhaps, memory_cols)

    colors = plt.cm.tab20.colors[:len(id_cols)+1]
    colors = colors[:5] + colors[6:]

    fig, axs = plt.subplots(2, 2, figsize=(20, 6), sharex=True)

    grouped_violin(axs[0][0], cenhaps, runtime_data, runtime_cols, colors)
    axs[0][0].set_ylabel("Runtime", fontsize=13)

    grouped_violin(axs[0][1], cenhaps, memory_data, memory_cols, colors)
    axs[0][1].set_ylabel("Memory", fontsize=13)

    grouped_violin(axs[1][0], cenhaps, id_data, id_cols, colors)
    axs[1][0].set_ylabel("Identity", fontsize=13)
    axs[1][0].set_xlabel("Truth cenhap", fontsize=13)

    grouped_violin(axs[1][1], cenhaps, cor_data, cor_cols, colors)
    axs[1][1].set_ylabel("Correctness", fontsize=13)
    axs[1][1].set_xlabel("Truth cenhap", fontsize=13)

    panel_letter(axs[0][0], "a")
    panel_letter(axs[0][1], "b")
    panel_letter(axs[1][0], "c")
    panel_letter(axs[1][1], "d")

    legend_from_columns(fig, id_cols, colors)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(outdir, name), dpi=300)
    plt.close(fig)

def sup_dists(rows, outdir):

    cenhaps = sorted({r["Truth cenhap"] for r in rows})
    rows_by_cenhap = {cenhap : [r for r in rows if r["Truth cenhap"] == cenhap]
                      for cenhap in cenhaps}
    cmap = plt.cm.tab10
    colors = {c: cmap(i) for i, c in enumerate(cenhaps)}

    plt.figure(figsize=(6, 6))

    for cenhap, cenhap_rows in rows_by_cenhap.items():
        plt.scatter([r["Minimum graph distance"] for r in cenhap_rows],
                    [r["Minimum sampled distance"] for r in cenhap_rows],
                    s=12, alpha=0.7, color=colors[cenhap], label=cenhap)

    plt.axline((0, 0), slope=1)

    plt.xlabel("Minimum graph distance", fontsize=13)
    plt.ylabel("Minimum sampled distance", fontsize=13)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sup_dists.png"), dpi=300)
    plt.close()

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

    sup_violin(rows, header, args.output_dir, True,
               "sup_real_02.png", lambda d: d < 0.2)

    sup_violin(rows, header, args.output_dir, True,
               "sup_real_02_04.png", lambda d: 0.2 < d < 0.4)

    sup_violin(rows, header, args.output_dir, True,
               "sup_real_04.png", lambda d: d > 0.4)
    
    sup_violin(rows, header, args.output_dir, False,
               "sup_sim_02.png", lambda d: d < 0.2)

    sup_violin(rows, header, args.output_dir, False,
               "sup_sim_02_04.png", lambda d: 0.2 < d < 0.4)

    sup_violin(rows, header, args.output_dir, False,
               "sup_sim_04.png", lambda d: d > 0.4)

    sup_dists(rows, args.output_dir)

    sup_correctness(rows, args.output_dir)