import os
import pandas as pd

def find_fna_files(root_dir):
    fna_files = {}
    for dirpath, _, filenames in os.walk(root_dir):
        for fn in filenames:
            if fn.endswith('.fna'):
                fna_files[fn] = os.path.join(dirpath, fn)
    return fna_files

def main(csv_file, root_dir):
    df = pd.read_csv(csv_file)
    filenames_set = set(df['filename'])  # Assuming CSV has a 'filename' column with base names

    all_fna_files = find_fna_files(root_dir)

    # Map only those files in CSV to their full paths
    df['filepath'] = df['filename'].map(all_fna_files)

    # Save updated CSV with paths appended
    output_csv = csv_file.replace('.csv', '_with_paths.csv')
    df.to_csv(output_csv, index=False)
    print(f"Updated CSV saved to: {output_csv}")

    # Example: Run GTDB-tk classify on filtered list (default parameters)
    selected_files = df.loc[df['filepath'].notnull(), 'filepath'].tolist()
    for genome_file in selected_files:
        # Construct and print command for GTDB-tk (customize as needed)
        cmd = f"gtdbtk classify_wf --genome_dir {os.path.dirname(genome_file)} --out_dir gtdbtk_output --cpus 4"
        print(f"Run GTDB-tk with:\n{cmd}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py <csv_file_with_filenames> <root_directory_to_search>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
