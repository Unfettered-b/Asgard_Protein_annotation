#!/usr/bin/env python3
# ================================================
# Filter Asgard Genomes Based on CheckM2 Results
# Completeness > 50%, Contamination < 10%
# ================================================

import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import shutil


def filter_genomes(file, level):
    """Filter genomes analyzed using CheckM2 based on completeness and contamination."""

    print(f"\n[DEBUG] Reading CheckM2 results from: {file}")
    quality = pd.read_csv(file, sep='\t')
    print(f"[DEBUG] Columns found: {list(quality.columns)}")
    print(f"[DEBUG] Number of genomes before filtering: {len(quality)}")

    # Try to normalize the genome name column
    if "Genome" in quality.columns:
        quality.rename(columns={"Genome": "Name"}, inplace=True)
    elif "Name" not in quality.columns:
        raise KeyError("[ERROR] Could not find 'Genome' or 'Name' column in CheckM2 file!")

    # Convert numeric columns
    quality['Completeness'] = quality['Completeness'].astype(float)
    quality['Contamination'] = quality['Contamination'].astype(float)

    # Filtering flags
    quality['Contaminated'] = quality['Contamination'] > 10
    quality['Incomplete'] = quality['Completeness'] < int(level)
    quality['Selected'] = (~quality['Contaminated']) & (~quality['Incomplete'])

    # Genome Quality Score
    quality['GQS'] = quality['Completeness'] - 5 * quality['Contamination']

    filtered = quality[quality['Selected']]
    print(f"[DEBUG] Number of genomes passing filter (Completeness ≥ {level}, Contamination ≤ 10): {len(filtered)}")

    # Plot bar chart summary
    counts = {
        'Selected': quality['Selected'].sum(),
        'Contaminated': quality['Contaminated'].sum(),
        'Incomplete': quality['Incomplete'].sum()
    }

    plt.figure(figsize=(8, 6))
    plt.bar(counts.keys(), counts.values(), color=['blue', 'red', 'green'])
    plt.ylabel('Count')
    plt.title('Genome Quality Categories')
    plt.show()

    return filtered


def get_fna_files(root_dir):
    """Get all .fna files from root_dir with their full paths and names (without extension)."""

    print(f"\n[DEBUG] Scanning for .fna files in: {root_dir}")
    fna_files, fna_files_wo_root = [], []
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.lower().endswith('.faa'):
                fna_files.append(os.path.join(root, file))
                fna_files_wo_root.append(os.path.splitext(file)[0])

    fna_df = pd.DataFrame({"paths": fna_files, "Name": fna_files_wo_root})
    print(f"[DEBUG] Found {len(fna_df)} .faa files")
    return fna_df


def move_files(path_df, destination_path):
    """Move selected genomes to a destination directory."""
    os.makedirs(destination_path, exist_ok=True)

    print(f"\n[DEBUG] Copying {len(path_df)} files to {destination_path}")
    copied = 0
    for path in path_df.paths:
        filename = os.path.basename(path)
        dest_path = os.path.join(destination_path, filename)
        try:
            shutil.copy2(path, dest_path)
            copied += 1
        except Exception as e:
            print(f"[WARNING] Could not move {path}: {e}")

    print(f"[DEBUG] Successfully copied {copied}/{len(path_df.paths)} files.")


def main(file):
    # Filter CheckM2 output
    level = 50
    filtered = filter_genomes(file, level)

    # Get .fna file paths
    root_path = "/home/anirudh/genomes/selected_genomes"
    path_df = get_fna_files(root_path)


    # Normalize names for merging
    filtered["Name"] = filtered["Name"].str.replace(".fna", "", regex=False).str.strip()
    path_df["Name"] = path_df["Name"].str.strip()

    # Merge filtered genomes with available paths
    print(f"\n[DEBUG] Merging filtered genomes with available .fna files...")
    merged = pd.merge(filtered, path_df, on="Name", how='inner')
    print(f"[DEBUG] Number of genomes matched between CheckM2 and .fna files: {len(merged)}")

    # Save merged list
    filtered_output = f'/home/anirudh/genomes/Asgard_genomes/Data/Filtered_genomes{level}.csv'
    merged.to_csv(filtered_output, index=False)
    print(f"[DEBUG] Filtered genome list exported to: {filtered_output}")

    # Move files
    final_dir = f'/home/anirudh/genomes/complete{level}/prokka'
    move_files(merged, final_dir)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python ./filter_genomes.py ./path/to/checkm2/output.tsv")
        print("Using default TSV file instead.")
        def_tsv = '/home/anirudh/genomes/Asgard_genomes/checkm2_asgard/quality_report.tsv'
    else:
        def_tsv = None

    file = sys.argv[1] if len(sys.argv) >= 2 else def_tsv
    main(file)
