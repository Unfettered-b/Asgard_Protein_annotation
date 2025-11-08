import pandas as pd
import os


inputdir = "/home/anirudh/genomes/asCOGs/results/"
clustersdir = os.path.join(inputdir, "selected")



countdb = pd.read_csv(os.path.join(clustersdir, "cluster_sizes.txt"), sep="\t", header=None, names=[ "size", "cluster"])
summary_file = os.path.join(inputdir, "asgard_clustering_summary.tsv")

countdb['size'] = pd.to_numeric(countdb['size'], errors='coerce')

print(countdb.head())


high_count_clusters = countdb[countdb['size'] >= 5]
num_clusters = high_count_clusters.shape[0]
num_proteins = high_count_clusters['size'].sum()

print(f"Number of clusters with >= 5 proteins: {num_clusters}")
print(f"Total number of proteins in these clusters: {num_proteins}")

with open(summary_file, "a") as f:
    f.write(f"Number of clusters with >= 5 proteins: {num_clusters}\n")
    f.write(f"Total number of proteins in these clusters: {num_proteins}\n")