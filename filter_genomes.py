# Program to remove genomes with completeness < 50 % and contamination > 10 %

import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import shutil


def filter_genomes(file):

    # Function to filter genomes analyzed using CheckM2 such that the completeness is >  50% and contamination is < 10%
    # file : output of CheckM2

    quality = pd.read_csv(file, sep='\t')

    print(f"Number of genomes: {len(quality)}")

    quality['Completeness'] = quality['Completeness'].astype(float)
    quality['Contamination'] = quality['Contamination'].astype(float)

    quality['Contaminated'] = quality['Contamination'] > 10
    quality['Incomplete'] = quality['Completeness'] < 50

    # Use bitwise OR for combining conditions element-wise
    quality['Selected'] = (~quality['Contaminated']) & (~quality['Incomplete'])


    quality['GQS'] = quality['Completeness'] - 5 * quality['Contamination']

    filtered = quality[quality['Selected']]


    # Plotting bar chart
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
    # Function to get all fna files from the root dir 
    # Creates a df with full path and the filename without extension to perform an inner merge with selected genomes df

    fna_files = []
    fna_files_wo_root = []
    for root, dirs, files in os.walk(root_dir):
        for file in files: 
            if file.lower().endswith('.fna'):
                fna_files.append(os.path.join(root, file))
                fna_files_wo_root.append(os.path.splitext(file)[0])

    fna_df = pd.DataFrame(
        {
            "paths" : fna_files,
            "Name": fna_files_wo_root
        }
    )
    return fna_df



def move_files(path_df, destination_path):
    # Function to move the selected genomes to one single directory for GTDB-TK classification
    # path_df : inner merged df on selected and the paths
    # destination_path : the selected genomes directory

    # creating the destination directory if it doesnt exist
    os.makedirs(destination_path, exist_ok= True)

    # Iterate through selected genomes and move the files
    for path in path_df.paths:
        filename = os.path.basename(path)
        dest_path = os.path.join(destination_path, filename)
        shutil.move(path, dest_path)

    print(f"Moved {len(path_df.paths)} files to {destination_path}")




def main(file):    
    # Filter outputs of CheckM2
    filtered = filter_genomes(file)
    

    # Get all .fna files from root_path
    root_path = "/mnt/d/genomes/Asgard_genomes"
    path_df = get_fna_files(root_path)

    merged = pd.merge(filtered, path_df, on="Name", how='inner')

    # print(merged)

    # Save the merged dataframe
    filtered_output = '/mnt/d/genomes/Asgard_genomes/Data/Filtered_genomes.csv'
    merged.to_csv(filtered_output, index=False)
    print(f'Filtered Genome list exported to {filtered_output}')

    # move the selected files to selected_genomes dir
    final_dir = '/mnt/d/genomes/selected_genomes'
    move_files(merged, final_dir)
    

    


    




if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Use python ./filter_genomes.py ./path/to/check2m/output.tsv")
        print('Using default tsv file')

        def_tsv = '/mnt/d/genomes/Asgard_genomes/checkm2_asgard/quality_report.tsv'

        
    file = sys.argv[1] if (len(sys.argv) >= 2) else def_tsv
    main(file)