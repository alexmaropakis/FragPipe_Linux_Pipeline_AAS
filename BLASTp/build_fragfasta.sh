#!/bin/bash
#SBATCH --job-name=mtp_pipeline
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/logs/mtp_pipeline_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/logs/mtp_pipeline_%j.err

set -euo pipefail

export PATH=$HOME/bin/ncbi-blast-2.17.0+/bin:$PATH
blastp -version

cd "$SLURM_SUBMIT_DIR"          # run from the dir holding the .py scripts
DEP=/scratch/maropakis.a/Dependencies
mkdir -p "$DEP/mtp_maps" "$DEP/FASTA_fragpipe"

# Stage 1: blast MTPs -> per-species CSVs (query id = full header + TMT tag)
python3 blast_mtps.py \
  --mtp-dir    "$DEP/FASTA_appended/" \
  --human-ref  "$DEP/FASTA/HUMAN.fasta" \
  --mouse-ref  "$DEP/FASTA/MOUSE_UP000000589_10090.fasta" \
  --out-dir    "$DEP/mtp_maps/" \
  --threads "${SLURM_CPUS_PER_TASK}"

# Stage 2: build FragPipe FASTAs from filtered CSVs
python3 build_fragpipe_fasta.py \
  --mtp-dir   "$DEP/FASTA_appended/" \
  --human-csv "$DEP/mtp_maps/human_filtered.csv" \
  --mouse-csv "$DEP/mtp_maps/mouse_filtered.csv" \
  --out-dir   "$DEP/FASTA_fragpipe/"

# Stage 3: sanity check — no duplicate headers in any FragPipe FASTA
fail=0
for f in "$DEP"/FASTA_fragpipe/*_fragpipe.fasta; do
  dups=$(grep '^>' "$f" | sort | uniq -d || true)
  if [ -n "$dups" ]; then
    echo "ERROR: duplicate headers in $f:" >&2
    echo "$dups" | head >&2
    fail=1
  fi
done
[ "$fail" -eq 0 ] && echo "All FragPipe FASTAs have unique headers." || exit 1
