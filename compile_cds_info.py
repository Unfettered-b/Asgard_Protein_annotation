import pandas as pd
import os
from collections import defaultdict as ddict
import sys
import shutil

# Takes input of the folder containing faa files having cds of all genomes filtered to a specific completion level 



def get_proteins(filepath):
    """
    Extract all predicted CDS from a .faa file.

    Parameters:
        filepath (str): Path to the .faa file.

    Returns:
        dict: Dictionary with keys 'header', 'cds_ids', and 'sequence'.
    """

    proteins = ddict(list)
    with open(filepath, "r") as f:
        header, cds_id, seq = None, None, ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    proteins['header'].append(header)
                    proteins['cds_ids'].append(cds_id)
                    proteins['sequence'].append(seq)  

                header = line[1:]
                cds_id = line[1:].split()[0]
                seq = ""
            else:
                seq += line

        if header is not None:
            proteins['header'].append(header)
            proteins['cds_ids'].append(cds_id)
            proteins['sequence'].append(seq) 
    
    return proteins



def build_protein_table(faabasepath, metadata_csv):
    # Read metadata table
    meta = pd.read_csv(metadata_csv, sep="\t")
    meta['full_name'] = meta['Assembly Accession'] + '_' + meta['Assembly Name']

    faa_files = os.listdir(faabasepath)

    records = []
    for faa in faa_files:
        print(f"Processing {faa} ...")
        proteins = get_proteins(os.path.join(faabasepath, faa))
        filename = os.path.basename(faa).split('_genomic')[0]

        # find organism name for this file
        organism_info = meta.loc[meta['full_name'] == filename]
        if organism_info.empty:
            org_name = None
            org_strain = None
            org_breed = None
            print(f"Couldnt find organism info for {filename}")
        else:
            org_name = organism_info.iloc[0]['Organism Name']
            org_strain = organism_info.iloc[0]['Organism Infraspecific Names Strain']
            org_breed = organism_info.iloc[0]['Organism Infraspecific Names Breed']

        for cds, header, seq in zip(proteins['cds_ids'], proteins['header'], proteins['sequence']):
            records.append({
                'organism_name': org_name,
                'breed': org_breed,
                'strain': org_strain,
                'cds_id': cds,
                'header': header,
                'genome_file': filename,              
                'sequence': seq,
            })
        
        print(f"found {len(records)} proteins")

    df = pd.DataFrame(records)
    return df


def main():
        
    # define the input dirs
    completeness = sys.argv[1]  if len(sys.argv)>1 else 50
    genome_cds_folder = f"/home/anirudh/genomes/complete{completeness}/"
    genome_cds_paths = f"{genome_cds_folder}/prokka/"
    genomes_data = "/home/anirudh/genomes/Asgard_genomes/ncbi_dataset/ncbi_dataset.tsv"
    scripts_data = "/home/anirudh/genomes/scripts/data/"

    output = build_protein_table(genome_cds_paths, genomes_data)

    outfilename = os.path.join(genome_cds_folder, f"Proteins_genomes_cp{completeness}.csv")
    scripts_file = os.path.join(scripts_data, f"Proteins_genomes_cp{completeness}.csv")
    output.to_csv(outfilename)
    shutil.copy2(outfilename, scripts_file)

    print(f"Full Proteins file for completeness {completeness} has been saved to {outfilename}")


if __name__ == "__main__":
    main()








