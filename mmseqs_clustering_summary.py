#!/usr/bin/env python3
import sys
import os
import pandas as pd

if len(sys.argv) < 2:
    print("Usage: python generate_summary.py <output_dir>")
    sys.exit(1)

outdir = sys.argv[1]
assigned_dir = os.path.join(outdir, "assigned")
unassigned_dir = os.path.join(outdir, "unassigned")
denovo_dir = os.path.join(outdir, "denovo")
searches_dir = os.path.join(outdir, "searches")

print("[8] Generating summary statistics...")

def count_lines(filepath, grep_fasta=False):
    """Count lines in file, optionally counting only FASTA headers."""
    if not os.path.isfile(filepath):
        return 0
    if grep_fasta:
        return sum(1 for line in open(filepath, "r", encoding="utf-8", errors='ignore') if line.startswith(">"))
    return sum(1 for _ in open(filepath, "r", encoding="utf-8", errors='ignore'))

def count_unique_headers(fasta_file):
    """Count unique headers in FASTA file."""
    if not os.path.isfile(fasta_file):
        return 0
    headers = set()
    with open(fasta_file, "r", encoding="utf-8", errors='ignore') as f:
        for line in f:
            if line.startswith(">"):
                headers.add(line.strip())
    return len(headers)

# --- Basic counts ---
# Total query proteins from original FASTA or query DB
query_fasta = os.path.join(searches_dir, "query_db")
if os.path.exists(query_fasta):
    # Try to get sequence count from MMseqs2 DB
    try:
        result = os.popen(f"mmseqs stats {query_fasta} 2>/dev/null | grep 'Entries:'").read()
        total_query = int(result.split(":")[1].strip()) if "Entries:" in result else 0
    except:
        total_query = 0
else:
    # Fallback: count from assigned and unassigned
    assigned_count = count_lines(os.path.join(assigned_dir, "assigned_ids.txt"))
    unassigned_count = count_lines(os.path.join(unassigned_dir, "unassigned_headers.txt"))
    total_query = assigned_count + unassigned_count

# Assigned proteins
assigned_count = count_lines(os.path.join(assigned_dir, "assigned_ids.txt"))

# Unassigned proteins
unassigned_headers_file = os.path.join(unassigned_dir, "unassigned_headers.txt")
if os.path.exists(unassigned_headers_file):
    unassigned_count = count_lines(unassigned_headers_file)
else:
    unassigned_count = count_lines(os.path.join(unassigned_dir, "unassigned.faa"), grep_fasta=True)

# Filtered unassigned
filtered_unassigned_count = count_lines(os.path.join(unassigned_dir, "unassigned_60aa.faa"), grep_fasta=True)

# --- Assigned clusters ---
assigned_tsv = os.path.join(searches_dir, "ref_search.tsv")
if os.path.isfile(assigned_tsv):
    try:
        df_assigned = pd.read_csv(assigned_tsv, sep="\t", header=None, usecols=[1])
        assigned_clusters = df_assigned[1].nunique()
    except:
        assigned_clusters = 0
else:
    assigned_clusters = 0

# --- De novo clusters ---
denovo_tsv = os.path.join(denovo_dir, "denovo_clu.tsv")
if os.path.isfile(denovo_tsv):
    try:
        df_denovo = pd.read_csv(denovo_tsv, sep="\t", header=None, names=["cluster", "protein"])
        denovo_clusters = df_denovo["cluster"].nunique()
        denovo_proteins_clustered = df_denovo["protein"].nunique()
    except:
        denovo_clusters = 0
        denovo_proteins_clustered = 0
else:
    denovo_clusters = 0
    denovo_proteins_clustered = 0

# --- Percentages ---
assigned_pct = (assigned_count / total_query * 100) if total_query > 0 else 0
unassigned_pct = (unassigned_count / total_query * 100) if total_query > 0 else 0

# --- Print formatted summary ---
print("\n========================================")
print("📊 Asgard Clustering Summary Statistics")
print("========================================")
print(f"{'Category':<50}{'Count':>15}{'%':>10}")
print("-" * 75)
print(f"{'Total query proteins':<50}{total_query:>15}{'-':>10}")
print(f"{'Proteins assigned to reference clusters':<50}{assigned_count:>15}{assigned_pct:>10.2f}")
print(f"{'Unique reference clusters hit (assigned)':<50}{assigned_clusters:>15}{'-':>10}")
print(f"{'Unassigned proteins (no reference hit)':<50}{unassigned_count:>15}{unassigned_pct:>10.2f}")
print(f"{'Filtered unassigned proteins (≥60 aa)':<50}{filtered_unassigned_count:>15}{'-':>10}")
print(f"{'De novo clusters formed (unassigned)':<50}{denovo_clusters:>15}{'-':>10}")
print(f"{'Proteins in de novo clusters':<50}{denovo_proteins_clustered:>15}{'-':>10}")
print("========================================\n")

# --- Export TSV ---
summary_tsv = os.path.join(outdir, "asgard_clustering_summary.tsv")
summary_data = [
    ["Total_query_proteins", total_query, "-"],
    ["Proteins_assigned_to_reference_clusters", assigned_count, f"{assigned_pct:.2f}"],
    ["Unique_reference_clusters_hit", assigned_clusters, "-"],
    ["Unassigned_proteins", unassigned_count, f"{unassigned_pct:.2f}"],
    ["Filtered_unassigned_proteins_60aa", filtered_unassigned_count, "-"],
    ["De_novo_clusters_formed", denovo_clusters, "-"],
    ["Proteins_in_de_novo_clusters", denovo_proteins_clustered, "-"],
]

pd.DataFrame(summary_data, columns=["Category", "Count", "Percent"]).to_csv(
    summary_tsv, sep="\t", index=False
)

print(f"[DEBUG] Summary table saved to: {summary_tsv}")
print(f"✅ Done! All results and statistics are available in: {outdir}")