#!/bin/bash 

# Set input directory with default value
in_dir="${1:-/home/anirudh/genomes/asCOGs/results}"
# Set output directory for selected clusters
out_dir="$in_dir/selected"

echo "Input directory: $in_dir"
echo "Output directory: $out_dir"

# Create output directory if it doesn't exist
mkdir -p $out_dir

# Path to de novo clustering results
cluster_tsv="$in_dir/denovo/denovo_clu.tsv"

# Check if input files exist before proceeding
if [[ ! -f "$cluster_tsv" ]]; then
    echo "Error: Cluster file not found: $cluster_tsv"
    exit 1
fi

# Output format: <count> <cluster_id>
echo "Counting cluster sizes..."
awk '{print $1}' "$cluster_tsv" | sort | uniq -c | awk '{print $1, $2}' > "$out_dir/cluster_sizes.txt"

# Filter clusters with 5 or more members and extract cluster IDs
echo "Filtering large clusters (≥5 members)..."
awk '$1 >= 5 {print $2}' "$out_dir/cluster_sizes.txt" > "$out_dir/large_clusters.txt"

# Check if we have large clusters before proceeding
if [[ ! -s "$out_dir/large_clusters.txt" ]]; then
    echo "Warning: No large clusters found with ≥5 members"
    touch "$out_dir/denovo_reps_large.faa"  # Create empty file
else
    # Extract representative sequences for large clusters only
    # Using seqkit to grep sequences based on cluster IDs from large_clusters.txt
    echo "Extracting representatives for large clusters..."
    seqkit grep -f <(awk '{print $1}' "$out_dir/large_clusters.txt") \
        "$in_dir/denovo/denovo_reps.faa" > "$out_dir/denovo_reps_large.faa"
fi

# Combine reference cluster representatives with large de novo cluster representatives
# This creates the final set of representative sequences for large clusters
echo "Combining all representatives..."
cat "$in_dir/assigned/all_cluster_representatives.faa" \
    "$out_dir/denovo_reps_large.faa" > \
    "$out_dir/final_reps_large.faa"

# Verify the final file was created
if [[ ! -s "$out_dir/final_reps_large.faa" ]]; then
    echo "Error: Final representatives file is empty or not created"
    exit 1
fi

# Create MMseqs2 database from the final representative sequences
# This database will be used as the enriched database for ColabFold
echo "Creating MMseqs2 database..."
mmseqs createdb "$out_dir/final_reps_large.faa" "$out_dir/asgard_enriched_db"

# Print the size of the final enriched database
final_count=$(grep -c '^>' "$out_dir/final_reps_large.faa")
echo "========================================"
echo "Final enriched database: $final_count sequences"
echo "Database created at: $out_dir/asgard_enriched_db"
echo "========================================"

# Additional statistics
if [[ -s "$out_dir/cluster_sizes.txt" ]]; then
    total_clusters=$(wc -l < "$out_dir/cluster_sizes.txt")
    large_clusters=$(wc -l < "$out_dir/large_clusters.txt")
    echo "Cluster statistics:"
    echo "  Total clusters: $total_clusters"
    echo "  Large clusters (≥5 members): $large_clusters"
    echo "  Percentage large: $(echo "scale=2; $large_clusters * 100 / $total_clusters" | bc)%"
fi