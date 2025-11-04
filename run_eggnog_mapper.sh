#!/bin/bash

# =======================================================
# Annotation of CDS detected by Prokka using EggNOG-mapper
# Supports both Docker and Conda modes
# Parallelized using GNU parallel
# =======================================================

mode="$1"  # 1 = Docker mode, else Conda mode

# --- CONFIGURATION ---
BASE="/home/anirudh"
DATA_DIR="$BASE/genomes/eggnog"              # Local path to EggNOG DB
#INPUT_DIR="$BASE/genomes/selected_genomes/prokka_results" # when running fresh
INPUT_DIR="$BASE/genomes/missing_faa/miss2" # when running incomplete runs
OUTPUT_DIR_WSL="$BASE/emapper_results"       # WSL-native output
FINAL_OUTPUT_DIR="$BASE/genomes/selected_genomes/emapper_results"
CPU=8                                         # CPUs per job
JOBS=2                                        # Parallel jobs

DOCKER_IMAGE="nanozoo/eggnog-mapper:2.1.13--c16a7d2"


# --- Logging setup ---
LOG_DIR="$BASE/genomes/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/emapper_run_$(date '+%Y-%m-%d_%H-%M-%S').log"

# Redirect all output (stdout and stderr) to log file
exec > >(tee -a "$LOG_FILE") 2>&1

echo "==============================================="
echo " EggNOG-mapper pipeline started: $(date)"
echo " Mode: $([ "$mode" -eq 1 ] && echo Docker || echo Conda)"
echo " Log file: $LOG_FILE"
echo "==============================================="

# --- Ensure output directories exist ---
mkdir -p "$OUTPUT_DIR_WSL" "$OUTPUT_DIR_WSL/tmp" "$FINAL_OUTPUT_DIR"

# =======================================================
# Function to run a single file in Docker
# =======================================================
run_emapper() {
    local fasta="$1"
    local relpath
    relpath=$(realpath --relative-to="$INPUT_DIR" "$fasta")
    local outfile
    outfile="$(basename "$fasta" .faa).emapper"

    echo "[START] $fasta"

    docker run --rm \
        -v "$INPUT_DIR":/input:ro \
        -v "$OUTPUT_DIR_WSL":/output \
        -v "$DATA_DIR":/data:ro \
        --cpuset-cpus="0-15" \
        $DOCKER_IMAGE \
        emapper.py \
        -i "/input/$relpath" \
        --output "/output/$outfile" \
        --data_dir /data \
        --cpu $CPU \
        --override && \
        echo "[DONE] $fasta" || echo "[ERROR] $fasta"
}

# =======================================================
# Function to run a single file in Conda
# =======================================================
run_emapper_conda() {
    local fasta="$1"
    local relpath
    relpath=$(realpath --relative-to="$INPUT_DIR" "$fasta")
    local outfile
    outfile="$(basename "$fasta" .faa).emapper"

    echo "[START] $fasta"

    emapper.py \
        -i "$INPUT_DIR/$relpath" \
        --output "$OUTPUT_DIR_WSL/$outfile" \
        --data_dir "$DATA_DIR" \
        --cpu "$CPU" \
        --temp_dir "$OUTPUT_DIR_WSL/tmp" \
        --override \
        -m diamond && \
        echo "[DONE] $fasta" || echo "[ERROR] $fasta"
}

# --- Export functions and variables for parallel ---
export -f run_emapper run_emapper_conda
export INPUT_DIR OUTPUT_DIR_WSL DATA_DIR CPU DOCKER_IMAGE

# =======================================================
# Run all genomes in parallel
# =======================================================
if [ "$mode" -eq 1 ]; then
    echo "Running in Docker mode..."
    find "$INPUT_DIR" -type f -name "*.faa" | parallel -j "$JOBS" run_emapper {}
else
    echo "Running in Conda mode..."
    
    # Activate Conda environment once
    conda init
    conda activate eggnog
    
    find "$INPUT_DIR" -type f -name "*.faa" | parallel -j "$JOBS" run_emapper_conda {}
fi

# =======================================================
# Move results to Windows-mounted folder
# =======================================================
echo "[INFO] Moving all results to $FINAL_OUTPUT_DIR ..."
cp -r "$OUTPUT_DIR_WSL"/* "$FINAL_OUTPUT_DIR"/
echo "[INFO] Done!"
