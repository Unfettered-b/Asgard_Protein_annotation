#!/bin/bash

# script to search for particular domains in the obtained cds

# Take input of the specific profile to be searched for 

REF_PROF=$1
REF_HMM=$2


# Defining constant paths
OUTFILE="/home/anirudh/genomes/hmm_profiles/searches"
FINAL_DIR="/home/anirudh/genomes/asCOGs/results/selected"
REPS_FASTA="$FINAL_DIR/final_reps_large.faa"
REPS_DB="$FINAL_DIR/asgard_enriched_db"
CLUSTER_TSV="/home/anirudh/genomes/asCOGs/results/denovo/denovo_clu.tsv"
ALL_CDS="/home/anirudh/genomes/combined_proteins.faa"
