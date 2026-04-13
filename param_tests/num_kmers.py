#!/usr/bin/env python3
"""Make Panel B of --absent-score supplementary figure.

Take a file with graph-unique k-mer counts, and plot
the distribution among haplotypes.

Lines assumed to be sorted from small to large, as follows:
    HG00738#1#HG00738.1#0: 282 / 100686 kmers present (density 0.003)

./param_tests/num_kmers.py \
    -i /private/groups/cgl/jlsiren/centrolign/chr4.txt \
    -o plot_outputs/num_kmers.png
"""

import argparse # Command-line argument parsing

import matplotlib.pyplot as plt # Basic plotting

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', required=True,
                        help='Input k-mer stats file')
    parser.add_argument('-o', '--output-file', required=True,
                        help='Graph file to write to')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    kmer_counts = []
    first_hap_name = None
    with open(args.input_file) as file:
        for line in file:
            if '#' in line:
                parts = line.strip().split()
                if first_hap_name is None:
                    first_hap_name = parts[0].split('#')[2]
                kmer_counts.append(int(parts[1]))

    # Plot dotplot
    plt.figure(figsize=(8, 6))
    plt.plot(range(len(kmer_counts)), kmer_counts, 
             marker='.', linestyle='None')
    plt.plot(0, kmer_counts[0], label=first_hap_name, 
             marker='o', linestyle='None')
    plt.xticks([])

    plt.xlabel('Haplotype')
    plt.ylabel('# graph-unique k-mers')
    plt.legend()

    plt.tight_layout()
    plt.savefig(args.output_file, dpi=150)
    plt.close()
