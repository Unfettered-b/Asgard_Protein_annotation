#!/bin/bash
set -euo pipefail

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
query_faa="$1"
asgard_db="$2"
outdir="$3"

echo "========================================"
echo "[DEBUG] Starting Asgard Clustering Pipeline"
echo "[DEBUG] Query FASTA: $query_faa"
echo "[DEBUG] Asgard DB: $asgard_db"
echo "[DEBUG] Output directory: $outdir"
echo "========================================"

# --- Directory Setup ---
mkdir -p "$outdir/tmp"
# echo "[DEBUG] Created output directory structure under: $outdir/tmp"

# # --- Step 1: Create MMseqs2 DB for query ---
# echo "[DEBUG] Checking if input FASTA exists..."
# if [[ ! -f "$query_faa" ]]; then
#   echo "[ERROR] Input FASTA file not found: $query_faa"
#   exit 1
# fi

# echo "[DEBUG] Running: mmseqs createdb \"$query_faa\" \"$outdir/query_db\""
# mmseqs createdb "$query_faa" "$outdir/query_db"
# echo "[DEBUG] Query DB created at: $outdir/query_db"

# # --- Step 2: Assign query proteins to reference Asgard clusters ---
# echo "[1] Assigning query proteins to Asgard reference clusters..."
# echo "[DEBUG] Command: mmseqs search \"$outdir/query_db\" \"$asgard_db\" \"$outdir/ref_search\" \"$outdir/tmp\" -e 0.001 -s 9 --cov-mode 0 -c 0.8"
# mmseqs search "$outdir/query_db" "$asgard_db" "$outdir/ref_search" "$outdir/tmp" -e 0.001 -s 9 --cov-mode 0 -c 0.8
# mmseqs convertalis "$outdir/query_db" "$asgard_db" "$outdir/ref_search" "$outdir/ref_search.tsv"
# echo "[DEBUG] Reference search results saved to: $outdir/ref_search.tsv"
# --- Step 3: Identify unassigned proteins ---
echo "[2] Extracting unassigned proteins..."

# 1. Extract IDs of sequences that have hits
mmseqs convertalis "$outdir/query_db" "$asgard_db" "$outdir/ref_search" "$outdir/assigned.tsv" --format-output query,target

# 2. Get unique query IDs (sequences that had hits)
cut -f1 "$outdir/assigned.tsv" | sort -u > "$outdir/assigned_ids.txt"

# 3. Extract all query sequence IDs
mmseqs convert2fasta "$outdir/query_db" "$outdir/all_query_tmp.faa"
grep "^>" "$outdir/all_query_tmp.faa" | sed 's/^>//' > "$outdir/all_query_ids.txt"
rm "$outdir/all_query_tmp.faa"

# 4. Get IDs that were not assigned
grep -Fxv -f "$outdir/assigned_ids.txt" "$outdir/all_query_ids.txt" > "$outdir/unassigned_ids.txt"

# 5. Create a sub-DB with unassigned sequences
mmseqs createsubdb "$outdir/unassigned_ids.txt" "$outdir/query_db" "$outdir/unassigned_proteins_db"

# 6. Convert to FASTA
mmseqs convert2fasta "$outdir/unassigned_proteins_db" "$outdir/unassigned.faa"

echo "[DEBUG] Unassigned proteins saved to: $outdir/unassigned.faa"

# --- Step 4: Filter unassigned proteins (≥60 aa) ---
echo "[3] Filtering unassigned proteins (≥60 aa)..."
awk '/^>/ {if(seq) print seq; print; seq=""} /^[^>]/ {seq=seq$0} END{if(seq) print seq}' "$outdir/unassigned.faa" \
| awk 'BEGIN{RS=">"; FS="\n"} length($2) >= 60 {print ">"$0}' > "$outdir/unassigned_60aa.faa"

echo "[DEBUG] Filtered unassigned sequences saved to: $outdir/unassigned_60aa.faa"
echo "[DEBUG] File size check:"
ls -lh "$outdir/unassigned_60aa.faa"

# --- Step 5: Cluster unassigned proteins de novo ---
echo "[4] Clustering unassigned proteins de novo..."
mmseqs createdb "$outdir/unassigned_60aa.faa" "$outdir/unassigned_60aa_db"
mmseqs cluster "$outdir/unassigned_60aa_db" "$outdir/denovo_clu" "$outdir/tmp" --min-seq-id 0.2 -c 0.5
mmseqs createtsv "$outdir/unassigned_60aa_db" "$outdir/unassigned_60aa_db" "$outdir/denovo_clu" "$outdir/denovo_clu.tsv"
echo "[DEBUG] De novo clustering complete."

# --- Step 6: Build sequence profiles ---
echo "[5] Building sequence profiles..."
mmseqs createseqfiledb "$outdir/unassigned_60aa_db" "$outdir/denovo_clu" "$outdir/denovo_clu_seq"
mmseqs result2profile "$outdir/unassigned_60aa_db" "$outdir/unassigned_60aa_db" "$outdir/denovo_clu" "$outdir/denovo_profiles"
echo "[DEBUG] Sequence profiles generated."

# --- Step 7: Get representative sequences ---
echo "[6] Selecting representative sequences..."
mmseqs search "$outdir/unassigned_60aa_db" "$outdir/denovo_profiles" "$outdir/cluster_vs_profile" "$outdir/tmp"
mmseqs convertalis "$outdir/unassigned_60aa_db" "$outdir/denovo_profiles" "$outdir/cluster_vs_profile" "$outdir/cluster_vs_profile.tsv"
echo "[DEBUG] Representative sequences identified."

# --- Step 8: Merge representatives ---
echo "[7] Merging all representative sequences..."
rep_pattern="${asgard_db}_rep.*.faa"
echo "[DEBUG] Looking for Asgard representative files matching: $rep_pattern"
if compgen -G "$rep_pattern" > /dev/null; then
  cat $rep_pattern "$outdir/unassigned_60aa.faa" > "$outdir/all_cluster_representatives.faa"
  echo "[DEBUG] Merged file created: $outdir/all_cluster_representatives.faa"
else
  echo "[WARNING] No Asgard representative FASTA files found matching: $rep_pattern"
  cp "$outdir/unassigned_60aa.faa" "$outdir/all_cluster_representatives.faa"
fi

echo "✅ Done! All results saved in: $outdir"
