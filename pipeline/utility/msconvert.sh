#!/usr/bin/env bash
#SBATCH --job-name=raw2mzml
#SBATCH --partition=short
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=04:00:00


python /home/$USER/scripts/search_gen/msconvert.py /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_1 --plex cortex_1 --spectra-root /scratch/$USER/spectra
