# Check the actual cluster size distribution
in_dir="/home/anirudh/genomes/asCOGs/results"
out_dir="$in_dir/selected"

echo "=== CLUSTER SIZE ANALYSIS ==="
echo "Total clusters: $(wc -l < $out_dir/cluster_sizes.txt)"
echo ""

# Check the actual distribution
echo "Cluster size distribution:"
awk '{print $1}' "$out_dir/cluster_sizes.txt" | sort -n | uniq -c | head -20
echo ""

# Check the largest clusters
echo "Largest clusters:"
sort -nr "$out_dir/cluster_sizes.txt" | head -10