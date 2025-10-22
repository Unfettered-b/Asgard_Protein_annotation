#!/bin/bash
# Run GTDB-Tk classify_wf in batches using Docker

Defdir="/mnt/d/genomes/selected_genomes"
BATCH_SIZE=10  # Adjust based on RAM usage

if [ -z "$1" ]; then
  workdir="$Defdir"
else
  workdir="$1"
fi

outdir="$workdir/classified"
mkdir -p "$outdir"

# Get list of all genome files (assuming .fna/.fa)
genomes=("$workdir"/*.fna)
total=${#genomes[@]}

echo "Total genomes: $total"
echo "Running GTDB-Tk in batches of $BATCH_SIZE"

for ((i=0; i<total; i+=BATCH_SIZE)); do
    batch=("${genomes[@]:i:BATCH_SIZE}")
    
    # Create temporary batch folder
    batchdir="$workdir/batch_$i"
    mkdir -p "$batchdir"
    for g in "${batch[@]}"; do
        cp "$g" "$batchdir/"
    done

    echo "Processing batch $((i/BATCH_SIZE+1)) with ${#batch[@]} genomes..."

    # Run Docker on the batch folder
    sudo docker run --rm \
      -v "$batchdir":/data \
      -v /mnt/d/genomes/Asgard_genomes/Data/release226:/refdata \
      --memory=64g --cpus=16 \
      ecogenomic/gtdbtk \
      classify_wf \
        --genome_dir /data \
        --out_dir /data/classified \
        --cpus 16 \

    # Move results to main output folder
    cp "$batchdir/classified/"* "$outdir/"

    # Clean up batch folder
    rm -rf "$batchdir"
done

echo "All batches processed. Final results in $outdir"
