#!/bin/bash 

NEW_GENOME_FOLDER=$1

# extract and verify

for genome in "$input_dir"/*.zip; do
    extractNverify.sh "$genome"

# run checkM2
    Genome_QA.sh "$input_dir/$genome"


# prokka
    


# add to the combined proteins list


# add to compile cds



# program to add new genomes to the exisiting db
