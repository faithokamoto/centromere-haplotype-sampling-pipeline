#!/usr/bin/env python3

import os
import csv
import argparse
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-tsv', help='TSV with data to plot', required=True)
    parser.add_argument('-o', '--output-dir', help='Directory to write to', required=True)
    parser.add_argument('-n', '--haplotype-name', help='Haplotype to call out', required=True)
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


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    rows, _ = read_tsv(args.input_tsv)

    hg_depths = []
    other_depths = []

    for row in rows:
        names = row.get("Sampled haplotype names")
        depths = row.get("Sampled haplotype depths")

        if not names or not depths:
            continue

        # Split comma-separated values
        name_list = [x.strip() for x in str(names).split(",")]
        depth_list = [x.strip() for x in str(depths).split(",")]

        if len(name_list) == 0 or len(depth_list) == 0:
            continue

        try:
            first_name = name_list[0]
            first_depth = float(depth_list[0])
        except:
            continue

        # Cap depth at 30
        capped_depth = min(first_depth, 30)

        if first_name == args.haplotype_name:
            hg_depths.append(capped_depth)
        else:
            other_depths.append(capped_depth)

    # Plot grouped histogram
    bins = [val / 2 for val in range(0, 61)]  # 0–30 inclusive

    plt.figure(figsize=(8, 6))
    plt.hist(
        [hg_depths, other_depths],
        bins=bins,
        label=[args.haplotype_name, "Other"],
        alpha=0.7,
        edgecolor='black'
    )

    plt.xlabel("First haplotype depth (capped at 30)")
    plt.ylabel("Count")
    plt.title("Grouped histogram of first sampled haplotype depth")
    plt.legend()

    out_path = os.path.join(args.output_dir, "first_haplotype_depth_hist.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    print(f"Wrote plot to: {out_path}")