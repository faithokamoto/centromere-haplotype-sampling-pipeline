#!/usr/bin/env python3
"""Plot how well we found neighbor haplotypes.

Plot 
"""

from typing import Dict, List, Tuple

import matplotlib.pyplot as plt

DIST_THRESHOLD = 0.15
"""Maximum distance for a 'neighbor' haplotype."""
LAST_SAMPLE_RUN = 'Sampling 8 haplotypes'
"""The first thing logged during the last sampling run."""

# I/O file ptahs
DISTANCE_MATRIX = 'input_data/chr12_r2_QC_v2_centrolign_pairwise_distance.csv'
"""Three-column CSV of haplotype pairs and their distances."""
SAMPLE_TABLE = 'input_data/test_samples.txt'
"""Four-column CSV of (path name, version, haplotype name, optimal N)."""
LOG_FILE = 'log/test_output.txt'
OUTPUT_FILE = 'plot_outputs/neighbor_finding.png'

Matrix = Dict[str, Dict[str, float]]
"""{hap1: {hap2: distance}} dictionary."""
SampleTable = Dict[str, Tuple[str, int]]
"""{haplo name : (path name, optimal N)} dictionary."""
SampledHaplo = List[Tuple[str, float]]
"""List of (haplotype, score) tuples."""

def read_distance_matrix(dist_path: str) -> Matrix:
    """Read a distance matrix presented as a CSV of pairs."""
    distance_matrix = {}
    with open(dist_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            if len(parts) == 3:
                hap1, hap2, dist = parts
                dist = float(dist)

                # Just in case, we need to set up these entries
                if hap1 not in distance_matrix:
                    distance_matrix[hap1] = {}
                if hap2 not in distance_matrix:
                    distance_matrix[hap2] = {}

                distance_matrix[hap1][hap2] = dist
                distance_matrix[hap2][hap1] = dist
    return distance_matrix

def read_sample_table(sample_table_path: str) -> SampleTable:
    """Read a sample table and return a list of haplotype names."""
    sample_table = dict()
    with open(sample_table_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            if len(parts) == 4:
                # path name, version, haplotype name, optimal N
                path_name, _, haplo_name, optimal_n = parts
                sample_table[haplo_name] = (path_name, int(optimal_n))
    return sample_table

def read_sampled_haplotypes(log_file: str, haplo_name: str) -> SampledHaplo:
    """Read haplotypes sampled for a given input in a log file."""
    sampled_haplotypes = []
    # How far have we read so far?
    found_haplo_name = False
    found_last_sample_run = False
    with open(log_file, 'r') as file:
        for line in file:
            # Start of this sample's log
            if 'Processing sample' in line and haplo_name in line:
                found_haplo_name = True
            # Start of this sample's last sampling run
            if found_haplo_name and LAST_SAMPLE_RUN in line:
                found_last_sample_run = True
            # A haplotype sampled in the last sampling run
            if found_last_sample_run and 'Selected haplotype' in line:
                # Selected haplotype <name> with score <score>
                haplo = line.split()[2]
                score = float(line.split()[5])
                sampled_haplotypes.append((haplo, score))
            # End of this sample's log
            if found_last_sample_run and 'autoindex' in line:
                break
    return sampled_haplotypes

def plot_score_drop(all_sampled_haplotypes: Dict[str, SampledHaplo],
                    all_optimal_ns: Dict[str, int], ax: plt.Axes) -> None:
    """Plot the score drop data for all haplotypes.

    For each haplotyp, normalize to the best score (y=0)
    and the optimal number to sample (x=0).
    Average Y values across the same X value.

    Parameters
    ----------
    all_sampled_haplotypes : Dict[str, SampledHaplo]
        A dictionary mapping haplotype names to their sampled haplotypes.
    all_optimal_ns : Dict[str, int]
        A dictionary mapping haplotype names to their optimal N values.
    ax : plt.Axes:
        The matplotlib Axes object to plot on.
    """

    # {x: [y1, y2, ...]}
    y_values = dict()
    for haplo_name, sampled_haplotypes in all_sampled_haplotypes.items():
        # Collect (x, y) data points
        score_drop_data = []

        # Reference points for normalization
        optimal_n = all_optimal_ns[haplo_name]
        best_score = sampled_haplotypes[0][1]

        for rank, (_, score) in enumerate(sampled_haplotypes):
            # Calculate normalized (x, y) points
            x = rank - optimal_n
            y = score - best_score
            # Remember points indexed by x value
            if x not in y_values:
                y_values[x] = []
            y_values[x].append(y)

    # Average Y values across the same X value
    avg_y_values = {x: sum(y_list) / len(y_list) 
                    for x, y_list in y_values.items()}
    # Plot a line from min X to max X
    sorted_x = sorted(y_values.keys())
    sorted_y = [avg_y_values[x] for x in sorted_x]
    ax.plot(sorted_x, sorted_y, marker='o', color='black')

    ax.set_xlabel('Dist from optimal N')
    ax.set_ylabel('Drop from max score')

def plot_nearest_dist(distance_matrix: Matrix,
                      all_sampled_haplotypes: Dict[str, SampledHaplo],
                      sample_table: Dict[str, SampledHaplo],
                      ax: plt.Axes) -> None:
    """Plot the dist to true neighbor vs. dist to first sampled.

    Parameters
    ----------
    distance_matrix : Matrix
        {hap1: {hap2: distance}} dictionary.
    all_sampled_haplotypes : Dict[str, SampledHaplo]
        A dictionary mapping haplotype names to their sampled haplotypes.
    sample_table : SampleTable
        A dictionary mapping haplotype names to their path name and optimal N.
    ax : plt.Axes:
        The matplotlib Axes object to plot on.
    """
    for haplo_name, sampled_haplotypes in all_sampled_haplotypes.items():
        path_name = sample_table[haplo_name][0]
        # Not all haplotypes appear in the distance matrix for some reason
        if path_name not in distance_matrix:
            continue
        
        # Find nearest distance matrix neighbor
        matrix_neighbor = (None, float('inf'))
        for neighbor, dist in distance_matrix[path_name].items():
            if dist < matrix_neighbor[1]:
                matrix_neighbor = (neighbor, dist)
        
        # Also look up distance to the first haplotype we sampled
        first_sampled = sampled_haplotypes[0][0]
        nearest_sampled_dist = distance_matrix[path_name][first_sampled]

        ax.plot(matrix_neighbor[1], nearest_sampled_dist, 
                marker='o', color='black')
    # Reference y=x line
    ax.axline((0, 0), color='red', slope=1)
    ax.set_xlabel('Dist to true neighbor')
    ax.set_ylabel('Dist to first sampled')

def plot_correctness(distance_matrix: Matrix,
                     all_sampled_haplotypes: Dict[str, SampledHaplo],
                     sample_table: SampleTable,
                     ax: plt.Axes) -> None:
    """Plot histogram of % correct neighbors found.

    Parameters
    ----------
    distance_matrix : Matrix
        {hap1: {hap2: distance}} dictionary.
    all_sampled_haplotypes : Dict[str, SampledHaplo]
        A dictionary mapping haplotype names to their sampled haplotypes.
    sample_table : SampleTable
        A dictionary mapping haplotype names to their path name and optimal N.
    ax : plt.Axes:
        The matplotlib Axes object to plot on.
    """
    hist_values = []
    for haplo_name, sampled_haplotypes in all_sampled_haplotypes.items():
        if sample_table[haplo_name][1] > 0:
            path_name = sample_table[haplo_name][0]
            optimal_n = sample_table[haplo_name][1]
            # If we have nothing to compare, we cn't be correct
            if optimal_n == 0 or not path_name in distance_matrix:
                continue

            # Count how many neighbors are in the matrix
            matrix_neighbors = set()
            for neighbor, dist in distance_matrix[path_name].items():
                if dist <= DIST_THRESHOLD:
                    matrix_neighbors.add(neighbor)

            # If there are no neighbors, there's nothing to be correct
            if not matrix_neighbors:
                continue

            # Count how many were sampled
            correct_neighbors = 0
            for sampled_haplo, _ in sampled_haplotypes[:optimal_n]:
                if sampled_haplo in matrix_neighbors:
                    correct_neighbors += 1

            hist_values.append(correct_neighbors / optimal_n * 100)

    ax.hist(hist_values, bins=range(0, 101, 10), edgecolor='black')
    ax.set_xlabel(f'% <{DIST_THRESHOLD} neighbors found')
    ax.set_ylabel('# other neighbors found')

def plot_neighbor_finding(distance_matrix: Matrix, 
                          sample_table: SampleTable,
                          log_file: str, output_file: str) -> None:
    """Plot how well we found neighbor haplotypes.

    Parameters
    ----------
    distance_matrix : Matrix
        {hap1: {hap2: distance}} dictionary.
    sample_table : SampleTable
        {haplo name : (path name, optimal N)} dictionary.
    log_file : str
        Path to the log file.
    output_file : str
        Path to save the plot.
    """

    # Read all necessary data
    all_sampled_haplotypes = {}
    all_optimal_ns = {}
    for haplo_name, (_, optimal_n) in sample_table.items():
        sampled_haplotypes = read_sampled_haplotypes(log_file, haplo_name)
        all_sampled_haplotypes[haplo_name] = sampled_haplotypes
        all_optimal_ns[haplo_name] = optimal_n

    # Make three-panel figure
    _, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))
    plot_score_drop(all_sampled_haplotypes, all_optimal_ns, axs[0][0])
    plot_nearest_dist(distance_matrix, all_sampled_haplotypes, 
                      sample_table, axs[0][1])
    plot_correctness(distance_matrix, all_sampled_haplotypes, sample_table, axs[1][0])
    axs[1][1].axis('off')
    plt.tight_layout()
    plt.savefig(output_file)

if __name__ == '__main__':
    distance_matrix = read_distance_matrix(DISTANCE_MATRIX)
    sample_table = read_sample_table(SAMPLE_TABLE)
    plot_neighbor_finding(distance_matrix, sample_table, LOG_FILE, OUTPUT_FILE)