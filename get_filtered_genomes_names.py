import pandas as pd
import matplotlib.pyplot as plt
import os

# --- Input paths ---
level = 70
tsv_path = "/home/anirudh/genomes/Asgard_genomes/ncbi_dataset/ncbi_dataset.tsv"
csv_path = f"/home/anirudh/genomes/Asgard_genomes/Data/Filtered_genomes{level}.csv"
output_path = f"/home/anirudh/genomes/Asgard_genomes/Data/Filtered_genomes{level}_names.csv"
plot_path = f"/home/anirudh/genomes/Asgard_genomes/Data/Organism_Distribution_{level}.png"

print(f"[DEBUG] Starting merge and analysis pipeline...")
print(f"[DEBUG] NCBI metadata TSV: {tsv_path}")
print(f"[DEBUG] Filtered genomes CSV: {csv_path}")

# --- Step 1: Read TSV (from NCBI datasets) ---
print(f"[DEBUG] Reading NCBI metadata...")
metadata = pd.read_csv(tsv_path, sep="\t", dtype=str)
print(f"[DEBUG] Columns in metadata: {list(metadata.columns)}")
print(f"[DEBUG] Total metadata entries: {len(metadata)}")

# --- Step 2: Create merged Name column ---
print(f"[DEBUG] Creating unified 'Name' column by merging 'Assembly Accession' and 'Assembly Name'...")
metadata["Name"] = (
    metadata["Assembly Accession"].str.strip()
    + "_"
    + metadata["Assembly Name"].str.replace(r"[^\w]+", "_", regex=True).str.strip()
)
print(f"[DEBUG] Example 'Name' values: {metadata['Name'].head(3).tolist()}")

# --- Step 3: Read the filtered genome CSV ---
print(f"[DEBUG] Reading filtered genomes CSV...")
filtered = pd.read_csv(csv_path, dtype=str)
filtered['Name'] = filtered['Name'].str.strip("_genomic")
print(f"[DEBUG] Columns in filtered CSV: {list(filtered.columns)}")
print(f"[DEBUG] Total filtered genomes: {len(filtered)}")

# --- Step 4: Merge on Name ---
print(f"[DEBUG] Merging filtered genomes with NCBI metadata on 'Name'...")
merged = pd.merge(filtered, metadata, on="Name", how="inner")
print(f"[DEBUG] Merge complete. Number of matched genomes: {len(merged)}")

# --- Step 5: Save output ---
os.makedirs(os.path.dirname(output_path), exist_ok=True)
merged.to_csv(output_path, index=False)
print(f"[DEBUG] ✅ Merged file saved to: {output_path}")

# --- Step 6: Organism Name statistics ---
if "Organism Name" in merged.columns:
    organism_counts = merged["Organism Name"].value_counts()
    print(f"\n📊 Organism Name Statistics:")
    print(f"[DEBUG] Total unique organisms: {organism_counts.shape[0]}")
    print(f"[DEBUG] Top 10 most represented organisms:\n{organism_counts.head(10)}")

    # --- Step 7: Plot top organisms ---
    top_n = 15
    print(f"[DEBUG] Plotting top {top_n} organisms by frequency...")
    plt.figure(figsize=(10, 6))
    organism_counts.head(top_n).plot(kind="barh")
    plt.xlabel("Number of Genomes")
    plt.ylabel("Organism Name")
    plt.title(f"Top {top_n} Organisms in Filtered Asgard Genomes")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    print(f"[DEBUG] 📈 Plot saved to: {plot_path}")
    plt.show()
else:
    print("[WARNING] Column 'Organism Name' not found in merged dataframe. Skipping plot.")
