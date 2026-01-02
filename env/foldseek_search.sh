#!/bin/bash
# Foldseek 2500 good PDBs vs uniprot50

INPUT_DIR="/home/anirudh/genomes/predicted"
TARGET_DB="/home/anirudh/genomes/databases/uniprot50"  # Fixed: databases/ (plural)
OUTPUT_DIR="/home/anirudh/genomes/foldseek_data"
FILTER_FILE="$INPUT_DIR/good_pdbs.txt"  # Paths to .pdb files

PDB_INPUT="$OUTPUT_DIR/foldseek_input"
RESULTS_DIR="$OUTPUT_DIR/foldseek_results"
TMP_DIR="$OUTPUT_DIR/foldseek_tmp"

# Create directories
mkdir -p "$PDB_INPUT" "$RESULTS_DIR" "$TMP_DIR"

# Clear old symlinks
rm -f "$PDB_INPUT"/*.pdb

# Create symlinks from good_pdbs.txt (one per line)
cd "$PDB_INPUT"
while IFS= read -r pdb_path; do
  if [[ -f "$INPUT_DIR/pdbs/$pdb_path.pdb" ]]; then
    ln -sf "$INPUT_DIR/pdbs/$pdb_path.pdb" "$(basename "$INPUT_DIR/pdbs/$pdb_path.pdb")"
  else
    echo "Warning: $INPUT_DIR/pdbs/$pdb_path.pdb not found"
  fi
done < "$FILTER_FILE"

# Verify
NUM_PDBS=$(ls "$PDB_INPUT"/*.pdb 2>/dev/null | wc -l)
echo "Created $NUM_PDBS symlinks"
du -sh "$PDB_INPUT"

if [[ $NUM_PDBS -eq 0 ]]; then
  echo "ERROR: No PDB symlinks created. Check good_pdbs.txt format."
  exit 1
fi

# RUN FOLDSEEK
echo "Running Foldseek on $NUM_PDBS PDBs vs uniprot50..."
foldseek easy-search \
  "$PDB_INPUT" \
  "$TARGET_DB" \
  "$RESULTS_DIR/result.tsv" \
  "$TMP_DIR" \
  --threads 32 \
  -s 9.5 \
  --sort-by-structure-bits 0 \
  --format-output query,target,evalue,tmscore,qalign,talign,alnlen,qid,tid,alnsize,mismatch,gapopen,raw

echo "Done! Results: $RESULTS_DIR/result.tsv"
echo "Significant hits (E-value < 0.001):"
awk -F'\t' '$3 < 0.001' "$RESULTS_DIR/result.tsv" | wc -l
head -5 "$RESULTS_DIR/result.tsv"
