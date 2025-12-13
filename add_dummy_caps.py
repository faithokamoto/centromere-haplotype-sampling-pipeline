#!/usr/bin/env python3
"""Add dummy start and end nodes to a GFA.

Input any number of GFA files. Node IDs will be changed,
both to avoid conflicts, and to add dummy nodes 1 and 2.
These start and end nodes will connect to all paths.

Input is a GFA. Haplotypes should be "P" paths.
Node IDs are assumed to be numeric.

Assumes GFA version 1.0: https://gfa-spec.github.io/GFA-spec/GFA1.html
"""

import argparse
from typing import Dict, List, Tuple

# Type aliases for clarity
Nodes = Dict[int, str]
"""node_id: sequence"""
Edge = Tuple[int, str, int, str, str]
"""(from_node, from_orient, to_node, to_orient, overlap)"""
Paths = Dict[str, List[Tuple[int, str]]]
"""path_name: [(node_id, orientation), ...]"""
Graph = Tuple[Nodes, List[Edge], Paths]
"""(nodes, [edge1, ...], paths)"""

# Magic numbers
DUMMY_START_ID = 1
DUMMY_END_ID = 2
HEADER_LINE = "H\tVN:Z:1.0"

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Add dummy start and end nodes to a GFA")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output GFA file with dummy nodes added")
    parser.add_argument("gfas", type=str,  nargs="+",
                         help="GFA files with haplotype paths")
    return parser.parse_args()

def read_gfa(gfa_file: str) -> Graph:
    """Read a GFA file and return a dictionary of haplotype paths.
    
    Parameters
    ----------
    gfa_file : str
        Path to the GFA file.

    Returns
    -------
    graph : Graph
        The graph as (nodes, edges, paths).
    """

    print("Reading GFA file:", gfa_file)
    nodes = dict()
    edges = []
    paths = dict()
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('H'):
                # ignore header lines for now
                pass
            # "S"/sequence lines are nodes
            elif line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = int(parts[1])
                seq = parts[2]
                nodes[node_id] = seq
            # "L"/link lines are edges
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                from_node = int(parts[1])
                from_orient = parts[2]
                to_node = int(parts[3])
                to_orient = parts[4]
                overlap = parts[5]
                edges.append((from_node, from_orient, to_node, to_orient, overlap))
            # "P"/path lines are paths
            elif line.startswith('P'):
                parts = line.strip().split('\t')
                paths[parts[1]] = [(int(node.strip('+-')), 
                                    '+' if node.endswith('+') else '-')
                                   for node in parts[2].split(',')]
            else:
                raise ValueError(f"Unrecognized GFA line: {line.strip()}")
    return nodes, edges, paths

def merge_gfas(gfa_files: List[str]) -> Graph:
    """Merge multiple GFA files into a single graph with unique node IDs.
    
    Parameters
    ----------
    gfa_files : List[str]
        List of paths to GFA files.

    Returns
    -------
    Graph
        The merged graph as (nodes, edges, paths).
    """

    merged_nodes = dict()
    merged_edges = []
    merged_paths = dict()
    current_max_ID = max(DUMMY_START_ID, DUMMY_END_ID)

    for gfa_file in gfa_files:
        nodes, edges, paths = read_gfa(gfa_file)
        id_mapping = dict()

        # Remap node IDs
        for old_id in nodes.keys():
            new_id = current_max_ID + 1
            id_mapping[old_id] = new_id
            merged_nodes[new_id] = nodes[old_id]
            current_max_ID += 1

        # Remap edges
        for edge in edges:
            from_node, from_orient, to_node, to_orient, overlap = edge
            new_edge = (id_mapping[from_node], from_orient,
                        id_mapping[to_node], to_orient, overlap)
            merged_edges.append(new_edge)

        # Remap paths
        for path_name, node_list in paths.items():
            new_node_list = [(id_mapping[node_id], orientation)
                             for node_id, orientation in node_list]
            merged_paths[path_name] = new_node_list

    return merged_nodes, merged_edges, merged_paths

def add_dummy_nodes(graph: Graph) -> Graph:
    """Add dummy start and end nodes to the graph.
    
    Parameters
    ----------
    graph : Graph
        The input graph as (nodes, edges, paths).

    Returns
    -------
    Graph
        The graph with dummy nodes added.
    """

    nodes, edges, paths = graph
    # Add dummy nodes
    nodes[DUMMY_START_ID] = 'N'  # Dummy start node
    nodes[DUMMY_END_ID] = 'N'  # Dummy end node

    # Connect dummy start node to all path starts
    for path_nodes in paths.values():
        first_node, first_orient = path_nodes[0]
        edges.append((DUMMY_START_ID, '+', first_node, first_orient, '*'))

    # Connect all path ends to dummy end node
    for path_nodes in paths.values():
        last_node, last_orient = path_nodes[-1]
        edges.append((last_node, last_orient, DUMMY_END_ID, '+', '*'))
    
    # Add dummy nodes into paths
    for path_name, path_nodes in paths.items():
        new_path = [(DUMMY_START_ID, "+")] + path_nodes + [(DUMMY_END_ID, "+")]
        paths[path_name] = new_path

    return nodes, edges, paths

def write_gfa(graph: Graph, output_file: str):
    """Write the graph to a GFA file.
    
    Parameters
    ----------
    graph : Graph
        The graph as (nodes, edges, paths).
    output_file : str
        Path to the output GFA file.
    """

    nodes, edges, paths = graph
    with open(output_file, 'w') as f:
        f.write(f"{HEADER_LINE}\n")
        for node_id, seq in sorted(nodes.items()):
            f.write(f"S\t{node_id}\t{seq}\n")
        for edge in edges:
            from_node, from_orient, to_node, to_orient, overlap = edge
            f.write(f"L\t{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t{overlap}\n")
        for path_name, node_list in paths.items():
            path_str = ','.join(f"{node_id}{orientation}" 
                                for node_id, orientation in node_list)
            f.write(f"P\t{path_name}\t{path_str}\t*\n")

if __name__ == "__main__":
    args = parse_args()
    print("Merging GFA files...")
    merged_graph = merge_gfas(args.gfas)
    print("Adding dummy nodes...")
    final_graph = add_dummy_nodes(merged_graph)
    print(f"Writing output to {args.output}...")
    write_gfa(final_graph, args.output)
