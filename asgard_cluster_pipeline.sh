#!/bin/bash
set +e -uo pipefail

# ================================
# Asgard Clustering Pipeline (Debug Mode)
# ================================
# Input:
#   1. Query proteins FASTA file
#   2. Reference MMseqs2 DB (asgard_db)
#   3. Output directory
# Usage:
#   ./asgard_cluster_pipeline.sh query_proteins.faa /path/to/asgard_db output_dir
# ================================

# --- Input Arguments ---
query_faa="${1:-/home/anirudh/genomes/selected_genomes/combined_proteins.faa}"
asgard_db="${2:-/home/anirudh/genomes/asCOGs/asgard_db}"
outdir="${3:-/home/anirudh/genomes/asCOGs/results}"

echo "========================================"
echo "[DEBUG] Starting Asgard Clustering Pipeline"
echo "[DEBUG] Query FASTA: $query_faa"
echo "[DEBUG] Asgard DB: $asgard_db"
echo "[DEBUG] Output directory: $outdir"
echo "========================================"

# --- Directory Setup ---
mkdir -p "$outdir/tmp"
mkdir -p "$outdir/assigned"
mkdir -p "$outdir/unassigned"
mkdir -p "$outdir/searches"
mkdir -p "$outdir/denovo"
echo "[DEBUG] Created output directory structure under: $outdir"

# # --- Step 1: Create MMseqs2 DB for query ---
# echo "[DEBUG] Checking if input FASTA exists..."
# if [[ ! -f "$query_faa" ]]; then
#   echo "[ERROR] Input FASTA file not found: $query_faa"
#   exit 1
# fi

# echo "[DEBUG] Running: mmseqs createdb \"$query_faa\" \"$outdir/searches/query_db\""
# mmseqs createdb "$query_faa" "$outdir/searches/query_db"
# echo "[DEBUG] Query DB created at: $outdir/searches/query_db"

# # --- Step 2: Assign query proteins to reference Asgard clusters ---
# echo "[1] Assigning query proteins to Asgard reference clusters..."
# echo "[DEBUG] Command: mmseqs search \"$outdir/searches/query_db\" \"$asgard_db\" \"$outdir/searches/ref_search\" \"$outdir/tmp\" -e 0.001 -s 9 --cov-mode 0 -c 0.8"
# mmseqs search "$outdir/searches/query_db" "$asgard_db" "$outdir/searches/ref_search" "$outdir/tmp" -e 0.001 -s 9 --cov-mode 0 -c 0.8

# # --- Step 3: Identify unassigned proteins using seqkit ---
# echo "[2] Extracting unassigned proteins..."

# Extract assigned IDs (first column only - query IDs)
mmseqs convertalis "$outdir/searches/query_db" "$asgard_db" "$outdir/searches/ref_search" "$outdir/assigned/assigned.tsv" --format-output query
cut -f1 "$outdir/assigned/assigned.tsv" | sort -u > "$outdir/assigned/assigned_ids.txt"

echo "[DEBUG] Assigned proteins: $(wc -l < "$outdir/assigned/assigned_ids.txt")"

# Use seqkit to extract unassigned sequences directly from original FASTA
echo "[DEBUG] Extracting unassigned sequences using seqkit..."
seqkit grep -v -f "$outdir/assigned/assigned_ids.txt" \
    "$query_faa" \
    > "$outdir/unassigned/unassigned.faa"

echo "[DEBUG] Unassigned proteins: $(grep -c '^>' "$outdir/unassigned/unassigned.faa")"

# Verify headers are unique
echo "[DEBUG] Verifying unique headers..."
unique_headers=$(grep '^>' "$outdir/unassigned/unassigned.faa" | sort -u | wc -l)
total_headers=$(grep -c '^>' "$outdir/unassigned/unassigned.faa")
if [[ "$unique_headers" -ne "$total_headers" ]]; then
    echo "[WARNING] Duplicate headers found! Creating unique headers..."
    awk 'BEGIN{counter=1} /^>/ {print ">unassigned_seq_" counter++; next} {print}' \
        "$outdir/unassigned/unassigned.faa" > "$outdir/unassigned/unassigned_fixed.faa"
    mv "$outdir/unassigned/unassigned_fixed.faa" "$outdir/unassigned/unassigned.faa"
    echo "[DEBUG] Fixed headers - now unique: $(grep -c '^>' "$outdir/unassigned/unassigned.faa")"
fi

# --- Step 4: Filter unassigned proteins (≥60 aa) ---
echo "[3] Filtering unassigned proteins (≥60 aa)..."
seqkit seq -m 60 --remove-gaps "$outdir/unassigned/unassigned.faa" > "$outdir/unassigned/unassigned_60aa.faa"

echo "[DEBUG] Filtered unassigned sequences: $(grep -c '^>' "$outdir/unassigned/unassigned_60aa.faa")"
echo "[DEBUG] Sample of filtered headers:"
# grep "^>" "$outdir/unassigned/unassigned_60aa.faa" | head -5

# --- Step 5: Cluster unassigned proteins de novo ---
echo "[4] Clustering unassigned proteins de novo..."
if [[ ! -s "$outdir/unassigned/unassigned_60aa.faa" ]]; then
    echo "[WARNING] No unassigned proteins to cluster. Creating empty files."
    touch "$outdir/denovo/denovo_clu.tsv"
    touch "$outdir/denovo/denovo_reps.faa"
else
    mmseqs createdb "$outdir/unassigned/unassigned_60aa.faa" "$outdir/unassigned/unassigned_60aa_db"
    
    # Cluster with reasonable parameters
    echo "[DEBUG] Running de novo clustering..."
    mmseqs cluster "$outdir/unassigned/unassigned_60aa_db" \
        "$outdir/denovo/denovo_clu" \
        "$outdir/tmp" \
        --min-seq-id 0.3 \
        -c 0.7 \
        --cov-mode 1
    
    mmseqs createtsv "$outdir/unassigned/unassigned_60aa_db" \
        "$outdir/unassigned/unassigned_60aa_db" \
        "$outdir/denovo/denovo_clu" \
        "$outdir/denovo/denovo_clu.tsv"

    mmseqs result2repseq "$outdir/unassigned/unassigned_60aa_db" \
        "$outdir/denovo/denovo_clu" \
        "$outdir/denovo/denovo_reps"

    mmseqs convert2fasta "$outdir/denovo/denovo_reps" \
        "$outdir/denovo/denovo_reps.faa"

    echo "[DEBUG] De novo clustering complete."
    num_denovo_clusters=$(grep -c '^>' "$outdir/denovo/denovo_reps.faa" 2>/dev/null || echo 0)
    num_input_seqs=$(grep -c '^>' "$outdir/unassigned/unassigned_60aa.faa")
    num_clustered_seqs=$(wc -l < "$outdir/denovo/denovo_clu.tsv" 2>/dev/null || echo 0)
    
    echo "[DEBUG] Input sequences: $num_input_seqs"
    echo "[DEBUG] Number of de novo clusters: $num_denovo_clusters"
    echo "[DEBUG] Sequences assigned to clusters: $num_clustered_seqs"
    
    # Check cluster distribution
    if [[ -f "$outdir/denovo/denovo_clu.tsv" ]]; then
        echo "[DEBUG] Cluster size distribution:"
        cut -f2 "$outdir/denovo/denovo_clu.tsv" | sort | uniq -c | sort -nr | head -5
    fi
fi

# --- Step 6: Merge representatives ---
echo "[5] Merging all representative sequences..."
# Look for Asgard representatives in the same directory as the database
asgard_dir=$(dirname "$asgard_db")
asgard_base=$(basename "$asgard_db")
rep_files=$(find "$asgard_dir" -name "${asgard_base}_rep*.faa" -o -name "${asgard_base}.rep*.faa" | head -1)

if [[ -n "$rep_files" && -s "$rep_files" ]]; then
  cat $rep_files > "$outdir/assigned/all_cluster_representatives.faa"
  echo "[DEBUG] Added Asgard representatives: $(grep -c '^>' $rep_files)"
else
  > "$outdir/assigned/all_cluster_representatives.faa"
  echo "[WARNING] No Asgard representative FASTA files found"
fi

# Add de novo representatives if they exist
if [[ -s "$outdir/denovo/denovo_reps.faa" ]]; then
  cat "$outdir/denovo/denovo_reps.faa" >> "$outdir/assigned/all_cluster_representatives.faa"
  echo "[DEBUG] Added de novo representatives: $(grep -c '^>' "$outdir/denovo/denovo_reps.faa")"
fi

echo "[DEBUG] Total representatives: $(grep -c '^>' "$outdir/assigned/all_cluster_representatives.faa" 2>/dev/null || echo 0)"

# Cleanup
rm -rf "$outdir/tmp"

echo "✅ Done! All results saved in: $outdir"

# The Python script runs automatically, no need to call it here
echo "[INFO] Pipeline completed. Genearting Summary."
python3 "/home/anirudh/genomes/scripts/mmseqs_clustering_summary.py" $outdir