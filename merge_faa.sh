#!/bin/bash

# --- Configuration ---
input_dir="/home/anirudh/genomes/selected_genomes/prokka_results"
output_csv="/home/anirudh/genomes/selected_genomes/all_faa_ids.csv"
output_fasta="/home/anirudh/genomes/selected_genomes/combined_proteins.faa"

# --- Initialize output files ---
echo "cds_id,source_file" > "$output_csv"
> "$output_fasta"   # truncate / create empty FASTA file

# --- Walk and process all .faa files ---
find "$input_dir" -type f -name "*.faa" | while read -r faa; do
  # Extract just the file name (e.g. GCA_001940725.faa)
  file_name="$(basename "$faa")"

  # Append sequence IDs + file info to CSV
  awk -v file="$file_name" '/^>/ {
      sub(/^>/, "", $0)
      print $1 "," file
  }' "$faa" >> "$output_csv"

  # Append all sequences to one combined FASTA
  cat "$faa" >> "$output_fasta"
done

echo "✅ Done! Outputs created:"
echo " - $output_csv"
echo " - $output_fasta"
