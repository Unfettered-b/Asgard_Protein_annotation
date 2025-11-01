# Python script to analyze the eggnog mapper results

import pandas as pd
import numpy as np
import sys
import os

mode = "debug"

def get_files(dir, extension):
    # finds all .annotation files in the given dir and its subdirectories
    annotation_files = [
        os.path.join(root, file)
        for root, dirs, files in os.walk(dir)
        for file in files
        if file.endswith(f"{extension}")
    ]

    print(f"Retrieved {extension} results; Found {len(annotation_files)} files")
    return annotation_files

def read_prokka_file(faa_file):
    """Read Prokka FAA file and extract gene IDs."""
    ids = []
    with open(faa_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])  # get header (first token)
    df = pd.DataFrame({"#query": ids})
    df["source_file"] = os.path.basename(faa_file)
    return df

def read_annotation_file(annotation_file):
    """Read EggNOG annotation file."""
    try:
        with open(annotation_file) as file:
                for i, line in enumerate(file):
                    if line.startswith("#query"):
                        header_line = i
                        break
                else:
                    raise ValueError("No header line starting with '#query' found")

            # Read the file starting from the header line
        df = pd.read_csv(annotation_file, sep="\t", header=header_line, low_memory=False)

        if "#query" not in df.columns:
            # Try to find the query column if renamed
            df.columns = [col.strip() for col in df.columns]
            possible = [c for c in df.columns if "query" in c.lower()]
            if possible:
                df = df.rename(columns={possible[0]: "#query"})
            else:
                raise ValueError("No #query column found.")
        df["source_file"] = os.path.basename(annotation_file)
        return df
    except Exception as e:
        print(f"Error reading {annotation_file}: {e}")
        return pd.DataFrame()

def compare_files(annotations, prokka_outputs):
    anno_bases = [os.path.basename(f).replace('.emapper.emapper.annotations', "") for f in annotations]
    faa_bases = [os.path.basename(f).replace('.faa', "") for f in prokka_outputs]

    common = set(anno_bases) & set(faa_bases)
    missing_in_ann = set(faa_bases) - set(anno_bases)

    print(len(common))
    return missing_in_ann
        
        
def merge_annotation_with_prokka(prokka_df, ann_df):
    """Merge Prokka and EggNOG data on #query."""
    if ann_df.empty:
        prokka_df["annotation_found"] = False
        return prokka_df.assign(**{col: "-" for col in ["COG_category", "Description"]})

    merged = pd.merge(prokka_df, ann_df, on="#query", how="left")
    for col in merged.select_dtypes(include="object").columns:
        merged[col] = merged[col].fillna("-")

    return merged


def process_pair(faa_file, ann_file):
    """Process one Prokka–EggNOG pair."""
    prokka_df = read_prokka_file(faa_file)
    ann_df = read_annotation_file(ann_file)
    return merge_annotation_with_prokka(prokka_df, ann_df)



def main(annotation_dir, prokka_dir):

    # Get all .annotation files from the results dir
    anns = get_files(annotation_dir, ".annotations")
    faas = get_files(prokka_dir, ".faa")

    missing = compare_files(anns, faas)
    with open("missing_faas.txt", "w") as out:
        out.write("\n".join(missing))

    combined = []
    for faa in faas:
        base = os.path.splitext(os.path.basename(faa))[0]
        # Try to find corresponding annotation file
        ann = next((f for f in anns if base in os.path.basename(f)), None)
        if ann:
            df = process_pair(faa, ann)
        else:
            df = read_prokka_file(faa)
            df["annotation_found"] = False
            df["COG_category"] = "-"
            df["Description"] = "-"
        combined.append(df)

    all_merged = pd.concat(combined, ignore_index=True)
    output = "/home/anirudh/genomes/prokka_eggnog_combined.csv"
    all_merged.to_csv(output, index=False)
    print(f"✅ Combined file saved: {output}")
    print(f"📊 Total entries: {len(all_merged)}")








if __name__=="__main__":
    defaultdir_annotation = "/home/anirudh/emapper_results"
    defaultdir_prokka = "/home/anirudh/genomes/selected_genomes/prokka_results"
    rootdir_annotation = sys.argv[1] if len(sys.argv) > 1 else defaultdir_annotation
    rootdir_prokka = sys.argv[2] if len(sys.argv) > 2 else defaultdir_prokka

    main(rootdir_annotation, rootdir_prokka)