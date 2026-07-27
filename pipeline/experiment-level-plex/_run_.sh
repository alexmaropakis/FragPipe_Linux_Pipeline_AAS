#!/usr/bin/env bash
#SBATCH --job-name=genworkflow
#SBATCH --partition=short
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=04:00:00

# If generating experiment-level runs, instead of running in terminal,
# list all the commands out here 
# This is where it helps to keep naming formats consistent

# Example
for tissue in cerebellum fat heart hippocampus kidney liver lung muscle spleen striatum; do
  python Alex_gen_fragpipe_experiments.py \
    --run        ${tissue}_Keele_2023 \
    --species    mouse \
    --workflow   /home/maropakis.a/scripts/search_gen/FragPipe/templates/TMT11_MS3_Val_TrypsinLysc.workflow \
    --sample-map /scratch/maropakis.a/Dependencies/sample_map/sample_map_keele_2023_${tissue}.xlsx \
    --experiment ${tissue}_1 /scratch/maropakis.a/MQ_raw/Keele_2023/${tissue}/b1 1 \
    --experiment ${tissue}_2 /scratch/maropakis.a/MQ_raw/Keele_2023/${tissue}/b2 2 \
    --fasta /scratch/maropakis.a/Dependencies/FASTA_database/saap_proteins_260727_mousedecoys.fasta \
    --out-dir    /scratch/maropakis.a/Frag_outputs/ \
    --spectra-root /scratch/maropakis.a/spectra
done