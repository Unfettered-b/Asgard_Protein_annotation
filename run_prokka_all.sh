#!/bin/bash

# ============================================
# Script: run_prokka_all_docker.sh
# Purpose: Annotate all .fna genomes in a folder using Prokka via Docker
# Author: Anirudh Bantwal Baliga
# ============================================

# Usage:
#   ./run_prokka_all_docker.sh /path/to/genomes /path/to/output
# Default output directory: ./prokka_results

# Exit on any error
set -e

# --- Input and Output Directories ---
input_dir=${1:-"./"}
output_dir=${2:-"/mnt/d/genomes/selected_genomes/prokka_results"}

# --- Create output directory if not exists ---
mkdir -p "$output_dir"

# --- Loop through all .fna files ---
for genome in "$input_dir"/*.fna; do
    # Skip if no .fna files exist
    [ -e "$genome" ] || { echo "No .fna files found in $input_dir"; exit 1; }

    # Extract base name (remove path and extension)
    base=$(basename "$genome" .fna)

    echo "Running Prokka (Docker) on: $base"

    # Dockerized Prokka run
    sudo docker run --rm \
        -v "$input_dir":/data/input \
        -v "$output_dir":/data/output \
        staphb/prokka:latest \
        prokka /data/input/"$base".fna \
        --outdir /data/output/"$base" \
        --prefix "$base" \
        --cpus 12 \
        --metagenome \
        --kingdom Archaea \
        --force

    echo "Finished annotation for: $base"
done

echo "All annotations completed. Results saved in: $output_dir"
