#!/usr/bin/env bash
#SBATCH --job-name=fragpipe
#SBATCH --partition=short
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/maropakis.a/Frag_outputs/logs/fp_%A_%a.out
#SBATCH --error=/scratch/maropakis.a/Frag_outputs/logs/fp_%A_%a.err
set -euo pipefail
export JAVA_HOME=$HOME/bin/jdk-17.0.18+8
export PATH=$JAVA_HOME/bin:$PATH

PLEX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/maropakis.a/Frag_outputs/plex_list.txt)

python3 /home/maropakis.a/scripts/FragPipe/run_plexes.py \
  --spectra-root /scratch/maropakis.a/spectra \
  --fasta-dir    /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
  --template-dir /home/maropakis.a/scripts/FragPipe/templates \
  --out-dir      /scratch/maropakis.a/Frag_outputs \
  --fragpipe-bin /home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe \
  --tools-folder /home/maropakis.a/fragpipe/fragpipe-24.0/tools \
  --only "$PLEX" --run

# sbatch --array=1-18 submit_fragpipe.sh to run
