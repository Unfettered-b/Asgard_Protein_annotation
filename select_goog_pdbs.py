# python script to create a txt file containting good pdbs

import pandas as pd
import sys



def main(file):
    data = pd.read_csv(file)
    print(f"Total number of sequences: {len(data)}")

    filtered_data = data[data["confidence"] >= 80]

    print(f"Number of sequences with confidence >= 80: {len(filtered_data)}")

    selected_ids = filtered_data["IDs"].tolist()
    with open("/home/anirudh/genomes/predicted/good_pdbs.txt", "w") as f:
        for id in selected_ids:
            f.write(f"{id}\n")
    print("Good pdb IDs written to good_pdbs.txt")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_good_pdbs_list.py <input_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    main(input_file)