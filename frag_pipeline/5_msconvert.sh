#!/usr/bin/env bash
#SBATCH --job-name=raw2mzml
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=$SPECTRA/logs/raw2mzml_%A_%a.log

set -uo pipefail

TRFP=~/thermoRawFileParser/ThermoRawFileParser
SPECTRA=/scratch/maropakis.a/spectra

cd TRFP 
mapfile -t PLEXES < <(find "$SPECTRA" -mindepth 1 -maxdepth 1 -type d ! -name logs | sort)
PLEX="${PLEXES[$((SLURM_ARRAY_TASK_ID - 1))]}"

for raw in "$PLEX"/*.raw; do
    [ -e "$raw" ] || continue
    base=$(basename "$raw" .raw)
    out="$PLEX/$base.mzML"
    [ -s "$out" ] && { echo "SKIP $base"; continue; }
    "$TRFP" -i="$raw" -o="$PLEX" -f=2 -l=3      # f=2 plain mzML; no -g
    [ -s "$out" ] && echo "OK $base" || echo "FAIL $base"
done

# Set `--array=1-N` to the number of plex folders (excluding `logs/`), then sbatch msconvert.sh
# Example: sbatch --array=1-18 msconvert.sh