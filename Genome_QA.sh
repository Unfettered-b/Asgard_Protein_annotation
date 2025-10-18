#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <path/to/extracted/folder>"
  exit 1
fi

input_dir="$1"
output_dir="$input_dir/checkm2_asgard"
# input_file="$input_dir/md5sum.txt"

mkdir -p "$output_dir"

# filtered_fna_list="fna_files_list.txt"
# > "$filtered_fna_list"


# Extract file paths, filter .fna files and save to list
# awk -v cwd="$input_dir" '{print cwd "/" $2}' "$input_file" | grep '\.fna$' > "$filtered_fna_list"


# Run CheckM2 predict on filtered list
checkm2 predict --threads 24 --input $(find $input_dir -type f -name "*.fna") --output-directory "$output_dir" --force


echo "CheckM2 run completed. Output at: $output_dir/quality_report.tsv"
