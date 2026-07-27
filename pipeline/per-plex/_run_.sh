#!/usr/bin/env bash
#SBATCH --job-name=buildplexes
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00

# If building multiple single-plexes at once, instead of running in terminal,
# write all commands to build the plexes here 

GEN=/home/maropakis.a/scripts/search_gen/Alex_gen_fragpipe.py
RAW=/scratch/maropakis.a/MQ_raw
SMAP=/scratch/maropakis.a/Dependencies/sample_map
FASTA=/scratch/maropakis.a/Dependencies/FASTA_fragpipe
TPL=/home/maropakis.a/scripts/search_gen/FragPipe/templates
OUT=/scratch/maropakis.a/Frag_outputs
SPECTRA=/scratch/maropakis.a/spectra

# gen <plex_token> <species> <workflow_file> <raw_subdir> <sample_map_file>
gen() {
  python3 "$GEN" "$RAW/$4" \
    --plex "$1" --species "$2" \
    --workflow "$TPL/$3" --sample-map "$SMAP/$5" \
    --fasta-dir "$FASTA" --out-dir "$OUT" --spectra-root "$SPECTRA"
}

# example
gen acgb1 human TMT10_MS3_Val.workflow Ping_2018/ACG/b1 sample_map_acgb1.xlsx

echo
echo "All plexes prepped. Submit them with:"
echo "  for s in $OUT/submit/submit_*.sh; do sbatch \"\$s\"; done"