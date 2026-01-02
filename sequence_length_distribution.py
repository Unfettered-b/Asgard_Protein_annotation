# program to plot the sequence length distribution

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys


def plot_seqlen_distribution(dataframe):
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(dataframe["length"], bins=30, color='skyblue', edgecolor='black')
    ax.set_title("Sequence Length Distribution")
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Frequency")
    plt.tight_layout()
    plt.show()


def data_wrangler(file_path):
    data = pd.read_csv(file_path)

    data["length"] = data["seqs"].apply(len)

    length800 = data[(data["length"] >= 800) & (data["length"] < 1500)]
    length1500 = data[data["length"] >= 1500]

    print(f"Number of sequences with length >= 800 and < 1500: {len(length800)}")
    print(f"Number of sequences with length >= 1500: {len(length1500)}")

    length_stats = data["length"].describe()
    print("Sequence Length Statistics:")
    print(length_stats)

    return data


def main(file):
    data = data_wrangler(file)

    plot_seqlen_distribution(data)



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sequence_length_distribution.py <path_to_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.isfile(input_file):
        print(f"Error: File '{input_file}' does not exist.")
        sys.exit(1)

    main(input_file)