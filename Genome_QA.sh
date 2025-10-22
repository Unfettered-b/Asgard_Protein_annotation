#!/bin/bash

# This program checks for the completeness and the contamination in the provided genomes 



if [ -z "$1" ]; then
  echo "Usage: $0 <path/to/extracted/folder>"
  exit 1
fi

input_dir="$1"
output_dir="$input_dir/checkm2_asgard"

mkdir -p "$output_dir"


# Run CheckM2 predict on filtered list
checkm2 predict --threads 24 --input $(find $input_dir -type f -name "*.fna") --output-directory "$output_dir" --force


echo "CheckM2 run completed. Output at: $output_dir/quality_report.tsv"
