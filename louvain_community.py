import networkx as nx
import csv
from community import community_louvain  # Install this library via `pip install python-louvain`
import argparse

def read_mash_results_to_graph(input_file):
    """
    Reads a Mash distance results file and creates a weighted graph.
    Assumes the file format: genome1, genome2, mash_distance, p_value, shared_hashes.
    """
    G = nx.Graph()

    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)  # Skip the header line
        for row in reader:
            genome1, genome2, mash_distance, _, _ = row
            mash_distance = float(mash_distance)

            # Calculate similarity as 1 - mash_distance (optional)
            similarity = 1 - mash_distance

            # Add edge to the graph
            G.add_edge(genome1, genome2, weight=similarity)

    return G

def detect_louvain_communities(graph, resolution):
    """
    Detects Louvain communities in a graph using the python-louvain library.
    Allows adjustment of the resolution parameter.
    """
    # Perform Louvain community detection with the specified resolution
    partition = community_louvain.best_partition(graph, weight='weight', resolution=resolution)
    return partition

def write_communities_to_file(partition, output_file):
    """
    Writes the Louvain community results to a file.
    """
    with open(output_file, 'w') as f:
        for genome, community in partition.items():
            f.write(f"{genome}\t{community}\n")
    print(f"Louvain community results written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Louvain communities from Mash distance results.")
    parser.add_argument("-i", "--input", required=True, help="Input Mash distance results file (tab-delimited).")
    parser.add_argument("-o", "--output", required=True, help="Output file for Louvain community results.")
    parser.add_argument("--resolution", type=float, default=1.0, help="Resolution parameter for Louvain algorithm (default: 1.0)")

    args = parser.parse_args()

    # Step 1: Read Mash results and create a graph
    graph = read_mash_results_to_graph(args.input)

    # Step 2: Detect Louvain communities with a specified resolution
    communities = detect_louvain_communities(graph, resolution=args.resolution)

    # Step 3: Write results to the output file
    write_communities_to_file(communities, args.output)
