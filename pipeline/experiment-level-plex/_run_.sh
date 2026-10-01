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

source activate  /projects/slavov/AM/envs/py39

python /home/$USER/scripts/Search_gen/FragPipe/gen_fragpipe_experiment_plex.py \
  --msconvert "$HOME/.conda/envs/pwiz/bin/msconvert" \
  --run FTLD \ # this is just what everything will be named prefix wise
  --workflow /home/$USER/scripts/Search_gen/FragPipe/templates/TMT-18-Val.workflow \
  --sample-map /scratch/$USER/Dependencies/sample_map/sample_map_ftld.xlsx \
  --experiment FTLD_1 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b1 1 \
  --experiment FTLD_2 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b2 2 \
  --experiment FTLD_3 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b3 3 \
  --experiment FTLD_4 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b4 4 \
  --experiment FTLD_5 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b5 5 \
  --experiment FTLD_6 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b6 6 \
  --experiment FTLD_7 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b7 7 \
  --experiment FTLD_8 /scratch/$USER/MQ_raw/Shrestha_2026/cohorts/mzXML/ftld/b8 8 \
  --fasta /scratch/$USER/Dependencies/FASTA_database/saap_proteins_260803_humandecoys.fasta \
  --out-dir /scratch/$USER/Frag_outputs \
  --spectra-root /scratch/$USER/spectra/
