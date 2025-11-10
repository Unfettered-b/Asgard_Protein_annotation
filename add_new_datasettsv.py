import pandas as pd

# Path to your existing CSV
file_path = "/home/anirudh/genomes/Asgard_genomes/ncbi_dataset/ncbi_dataset.tsv"

# Load the existing data
df = pd.read_csv(file_path, sep ="\t")

# Define the new row as a dictionary
new_row = {
    "Assembly Accession": "GCA_002505645.1",
    "Assembly Name": "ASM250564v1",
    "Organism Name": "Yibarchaeum umbracryptum",
    "Organism Infraspecific Names Breed": "updated recently",
    "Organism Infraspecific Names Strain": "updated recently",
    "Organism Infraspecific Names Cultivar": "updated recently",
    "Organism Infraspecific Names Ecotype": "updated recently",
    "Organism Infraspecific Names Isolate": "updated recently",
    "Organism Infraspecific Names Sex": "updated recently",
    "Annotation Name": "updated recently",
    "Assembly Level": "updated recently",
    "Assembly Release Date": "updated recently",
    "WGS project accession": "updated recently",
    "Assembly Stats Number of Scaffolds": "updated recently"
}

# Append and save back to file
df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
df.to_csv(file_path, sep='\t', index=False)

print("✅ Row appended successfully.")
