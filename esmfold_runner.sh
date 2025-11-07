#!/bin/bash
# ESMfold Docker Run Script
# Usage: ./esmfold_runner.sh [simple|custom]
# Not used since using ipynb

set -e

MODE=${1:-simple}

echo "========================================"
echo " ESMfold Structure Prediction"
echo "========================================"

# Configuration
IMAGE="docker.io/vivkr/esmfold:latest"
INPUT_DIR="/home/anirudh/genomes/asCOGs/results/selected"
OUTPUT_DIR="/home/anirudh/genomes/predicted"
FASTA_FILE="denovo_reps_large.faa"

mkdir -p "${OUTPUT_DIR}"

echo ""
echo "Configuration:"
echo "  Mode: ${MODE}"
echo "  Input: ${INPUT_DIR}/${FASTA_FILE}"
echo "  Output: ${OUTPUT_DIR}"
echo ""

if [ "${MODE}" == "simple" ]; then
    echo "Running ESMfold in SIMPLE mode..."
    echo ""

    docker run --gpus all --rm \
        -v "${INPUT_DIR}:/input:ro" \
        -v "${OUTPUT_DIR}:/output" \
        "${IMAGE}" \
        esmfold_infer.py \
        --fasta "/input/${FASTA_FILE}" \
        --outdir "/output"

elif [ "${MODE}" == "custom" ]; then
    echo "Running ESMfold in CUSTOM mode (Python script)..."
    echo ""

    docker run --gpus all --rm \
        -v "${INPUT_DIR}:/data/input:ro" \
        -v "${OUTPUT_DIR}:/data/output" \
        -v "$(pwd)/scripts:/scripts:ro" \
        "${IMAGE}" \
        python3 /scripts/esmfold_inference.py

else
    echo "Invalid mode: ${MODE}"
    echo "Usage: ./esmfold_runner.sh [simple|custom]"
    exit 1
fi

echo ""
echo "========================================"
echo "✓ ESMfold prediction completed!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "========================================"
