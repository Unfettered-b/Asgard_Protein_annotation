# Program to remove genomes with completeness < 50 % and contamination > 10 %

import pandas as pd
import matplotlib.pyplot as plt
import sys

def main(file):
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

    filtered_output = '/mnt/d/genomes/Asgard_genomes/Data/Filtered_genomes.csv'
    filtered.to_csv(filtered_output, index=False)
    print(f'Filtered Genome list exported to {filtered_output}')

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




if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Use python ./filter_genomes.py ./path/to/check2m/output.tsv")
        print('Using default tsv file')

        def_tsv = '/mnt/d/genomes/Asgard_genomes/checkm2_asgard/quality_report.tsv'

        
    file = sys.argv[1] if (len(sys.argv) >= 2) else def_tsv
    main(file)