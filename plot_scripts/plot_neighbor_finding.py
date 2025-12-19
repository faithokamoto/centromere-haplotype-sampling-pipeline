#!/usr/bin/env python3
"""Plot how well we found neighbor haplotypes.

Plot 
"""

from typing import Dict, List, Tuple

import matplotlib.pyplot as plt

DISTANCE_MATRIX = '/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/' \
                  'HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.txt'
DIST_THRESHOLD = 0.15
"""Maximum distance for a 'neighbor' haplotype."""
LAST_SAMPLE_RUN = 8

LOG_FILE = 'test_output.txt'
SAMPLE_TABLE = 'test_samples.txt'
OUTPUT_FILE = 'plot_outputs/neighbor_finding.png'

def read_distance_matrix(dist_path: str) -> Dict[str, Dict[str, float]]:
    """Read a distance matrix presented as a CSV of pairs.
    
    Parameters
    ----------
    dist_path : str
        Path to the distance matrix file.

    Returns
    -------
    Dict[str, Dict[str, float]]
        {hap1: {hap2: distance}} dictionary.
    """
    distance_matrix = {}
    with open(dist_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            if len(parts) == 3:
                hap1, hap2, dist = parts
                dist = float(dist)
                if hap1 not in distance_matrix:
                    distance_matrix[hap1] = {}
                if hap2 not in distance_matrix:
                    distance_matrix[hap2] = {}
                distance_matrix[hap1][hap2] = dist
                distance_matrix[hap2][hap1] = dist
    return distance_matrix

def find_neighbors(distance_matrix: Dict[str, Dict[str, float]], 
                   haplo_name: str, threshold: float) -> List[Tuple[str, float]]:
    """Find neighbors for a specific haplotype based on the distance threshold.

    Parameters
    ----------
    distance_matrix : Dict[str, Dict[str, float]]
        {hap1: {hap2: distance}} dictionary.
    haplo_name : str
        The name of the haplotype to find neighbors for.
    threshold : float
        The maximum distance to consider a haplotype a neighbor.

    Returns
    -------
    List[Tuple[str, float]]
        A list of (haplotype, distance) tuples for neighbors.
    """
    neighbors = []
    if haplo_name in distance_matrix:
        for neighbor, dist in distance_matrix[haplo_name].items():
            if dist <= threshold:
                neighbors.append((neighbor, dist))
    return neighbors

def read_sample_table(sample_table_path: str) -> Dict[str, Tuple[str, int]]:
    """Read a sample table and return a list of haplotype names.
    
    Parameters
    ----------
    sample_table_path : str
        Path to the sample table file.

    Returns
    -------
    Dict[str, Tuple[str, int]]
        {path name : (haplotype name, optimal N)} dictionary.
    """
    haplo_names = dict()
    with open(sample_table_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            if len(parts) == 4:
                path_name, _, haplo_name, optimal_n = parts
                haplo_names[path_name] = (haplo_name, int(optimal_n))
    return haplo_names

def read_sampled_haplotypes(log_file: str, 
                            haplo_name: str) -> List[Tuple[str, float]]:
    """Read the sampled haplotypes from a log file.
    
    Parameters
    ----------
    log_file : str
        Path to the log file.
    haplo_name : str
        The name of the haplotype to look for in the log file.
    
    Returns
    -------
    List[Tuple[str, float]]
        A tuple containing the list of sampled haplotypes and their scores.
    """
    sampled_haplotypes = []
    found_haplo_name = False
    found_last_sample_run = False
    with open(log_file, 'r') as file:
        for line in file:
            # Start of this sample's log
            if 'Processing sample' in line and haplo_name in line:
                found_haplo_name = True
            # Start of this sample's last sampling run
            if found_haplo_name and f'Sampling {LAST_SAMPLE_RUN} haplotypes' in line:
                found_last_sample_run = True
            # A haplotype sampled in the last sampling run
            if found_last_sample_run and 'Selected haplotype' in line:
                haplo = line.split()[2]
                score = float(line.split()[5])
                sampled_haplotypes.append((haplo, score))
            # End of this sample's log
            if found_last_sample_run and 'autoindex' in line:
                break
    return sampled_haplotypes
            
def get_score_drop_data(sampled_haplotypes: List[Tuple[str, float]], 
                        optimal_n: int) -> List[Tuple[float, float]]:
    """Get the score drop data for plotting.
    
    Normalize to the best score (y=0) and the optimal number (x=0)

    Parameters
    ----------
    sampled_haplotypes : List[Tuple[str, float]]
        A list of (haplotype, score) tuples for sampled haplotypes.
        Should be sorted already by score in descending order.
    optimal_n : int
        The optimal number of haplotypes.

    Returns
    -------
    List[Tuple[float, float]]
        A list of (normalized rank, normalized score drop) tuples.
    """
    best_score = sampled_haplotypes[0][1]
    score_drop_data = []
    for rank, (_, score) in enumerate(sampled_haplotypes):
        normalized_rank = rank - optimal_n
        normalized_score_drop = score - best_score
        score_drop_data.append((normalized_rank, normalized_score_drop))
    return score_drop_data

def plot_score_drop(all_sampled_haplotypes: Dict[str, List[Tuple[str, float]]],
                    all_optimal_ns: Dict[str, int], ax: plt.Axes) -> None:
    """Plot the score drop data for all haplotypes.

    Average Y values across the same X value.

    Parameters
    ----------
    all_sampled_haplotypes : Dict[str, List[Tuple[str, float]]]
        A dictionary mapping haplotype names to their sampled haplotypes.
    all_optimal_ns : Dict[str, int]
        A dictionary mapping haplotype names to their optimal N values.
    ax : plt.Axes:
        The matplotlib Axes object to plot on.
    """

    # {x: [y1, y2, ...]}
    y_values = dict()
    for haplo_name, sampled_haplotypes in all_sampled_haplotypes.items():
        optimal_n = all_optimal_ns[haplo_name]
        score_drop_data = get_score_drop_data(sampled_haplotypes, optimal_n)
        for x, y in score_drop_data:
            if x not in y_values:
                y_values[x] = []
            y_values[x].append(y)
    # Average Y values across the same X value
    avg_y_values = {x: sum(y_list) / len(y_list) 
                    for x, y_list in y_values.items()}
    sorted_x = sorted(y_values.keys())
    sorted_y = [avg_y_values[x] for x in sorted_x]
    # Plot the average Y values
    ax.plot(sorted_x, sorted_y, marker='o')
    ax.set_xlabel('Dist from optimal N')
    ax.set_ylabel('Drop from max score')

def plot_neighbor_finding(distance_matrix: Dict[str, Dict[str, float]], 
                          sample_table: Dict[str, Tuple[str, int]],
                          log_file: str, output_file: str) -> None:
    """Plot how well we found neighbor haplotypes.

    Parameters
    ----------
    distance_matrix : Dict[str, Dict[str, float]]
        {hap1: {hap2: distance}} dictionary.
    log_file : str
        Path to the log file.
    output_file : str
        Path to save the plot.
    """
    _, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))
    all_sampled_haplotypes = {}
    all_optimal_ns = {}
    for (haplo_name, optimal_n) in sample_table.values():
        sampled_haplotypes = read_sampled_haplotypes(log_file, haplo_name)
        all_sampled_haplotypes[haplo_name] = sampled_haplotypes
        all_optimal_ns[haplo_name] = optimal_n
    plot_score_drop(all_sampled_haplotypes, all_optimal_ns, axs[0][0])
    plt.savefig(output_file)

if __name__ == '__main__':
    distance_matrix = read_distance_matrix(DISTANCE_MATRIX)
    sample_table = read_sample_table(SAMPLE_TABLE)
    plot_neighbor_finding(distance_matrix, sample_table, LOG_FILE, OUTPUT_FILE)