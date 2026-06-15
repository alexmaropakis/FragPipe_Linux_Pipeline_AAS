# FragPipe in Linux Pipeline: a pipeline for running FragPipe in Linux (particularly for the Amino Acid Substitution Project)

[FragPipe](https://fragpipe.nesvilab.org/) is a comprehensive computational platform designed to analyze mass spectrometry-based proteomics data. It is powered by [MSFragger](https://msfragger.nesvilab.org/), an ultrafast proteomic search engine for both closed and open (wide precursor mass tolerance) peptide identification. It includes Percolator and the [Philosopher](https://nesvilab.github.io/philosopher/) toolkit for downstream post-processing of MSFragger search results, FDR filtering, and multi-experiment summary report generation. It also includes the [MSBooster](https://www.nature.com/articles/s41467-023-40129-9) module for deep-learning based rescoring of peptide identifications, and [Crystal-C](https://www.nesvilab.org/Crystal-C/) and [PTM-shepherd](https://github.com/Nesvilab/PTM-Shepherd) to aid interpretation of results from "open" and "mass offset" searches for PTMs.

This guide is targeted toward beginners planning to run FragPipe v24.0. It is optimized for [Northeastern University's Research Computing Cluster](https://rc.northeastern.edu/), which is hosted at the [Massachusetts Green High-Performance Computing Center (MGHPCC)](https://mghpcc.org/) in Holyoke, MA. It is particularly written to follow the AAS Pipeline defined in [Tsour et al., 2026](https://pubmed.ncbi.nlm.nih.gov/39253435/) and outlines a pipeline to generate FASTA files with custom peptides appended. 
* If not intending to use these custom FASTA files, skip step 1 and 2 and just generate the FASTA file + decoys in the FragPipe GUI.
* This pipeline is based off work done by Andrew Leduc, PhD at the Slavov Lab, located in [this repository](https://github.com/Andrew-Leduc/AAS_Evo), but optimized for my particular work on the AAS Project in the Slavov Lab. 

**Last updated:** 06-15-2026

## Overview
This repository automates preparation and execution of FragPipe searches, including:
- Construction of FragPipe-compatible FASTA databases
- Generation of annotation files from TMT sample maps
- Organization of spectra into per-plex directories
- RAW → mzML conversion using ThermoRawFileParser
- Headless FragPipe execution on HPC clusters
The pipeline also supports FASTA databases containing MTP entries generated from appended `*_MTP.fasta` files.

## Pipeline
```text
Reference FASTA
      ↓
1_PrepFASTA.py
      ↓
2_buildFragFASTA.py
      ↓
3_gen_annotations.py
      ↓
4_stage_spectra.py
      ↓
5_msconvert.sh
      ↓
7_run_plexes.py
      ↓
FragPipe
```

## Repository Structure

```text
FragPipe_Linux_Tutorial/

├── pipeline/
│   ├── 1_PrepFASTA.py
│   ├── 2_buildFragFASTA.py
│   ├── 3_gen_annotations.py
│   ├── 4_stage_spectra.py
│   ├── 5_msconvert.sh
│   ├── 6_submit_prep.txt
│   ├── 7_run_plexes.py
│   ├── _build_fragFASTA.sh
│   └── _submit_fragpipe.sh
│
├── workflow_templates/
│   ├── TMT10_MS3_Val.workflow
│   └── TMT16_Val.workflow
│
└── README.md
```

## Setting up FragPipe in Linux

The most recent FragPipe version (and older versions) can be accessed via an academic license [here](https://github.com/Nesvilab/FragPipe/releases). The official tutorial for using FragPipe can be found [here](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html).

Download the wanted zip file and put it in your home directory. The version of FragPipe this tutorial is optimized for is bundled in this repository with all dependencies. 

**Dependencies for FragPipe v24.0:**
* Require Java 11+
* Require MSFragger 4.4+
* Require IonQuant 1.11.18+
* Require diaTracer 2.2.1+
* Require Python 3.9, 3.10, or 3.11
* Upgrade Crystal-C to 1.5.10
* Upgrade MSBooster to 1.4.14
* Upgrade Philosopher to 5.1.3-RC9
* Upgrade PTM-Shepherd to 3.0.11
* Upgrade TMT-Integrator to 6.1.3
* Upgrade FragPipe-PDV to 1.5.6
* Upgrade FragPipe-SpecLib to 0.1.58

**Download Java 11+ runtime:**
```
cd ~/bin
wget https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.18%2B8/OpenJDK17U-jdk_x64_linux_hotspot_17.0.18_8.tar.gz
tar -xzf OpenJDK17U-jdk_x64_linux_hotspot_17.0.18_8.tar.gz
ls ~/java/jdk-17.0.18+8/bin
rm OpenJDK17U-jdk_x64_linux_hotspot_17.0.18_8.tar.gz

export JAVA_HOME=/home/maropakis.a/bin/jdk-17.0.18+8
export PATH=$JAVA_HOME/bin:$PATH
source ~/.bashrc

java -version
```

**Download FragPipe v24.0:**
```
mkdir /home/maropakis.a/fragpipe
cd fragpipe
unzip FragPipe-24.0-linux.zip
rm FragPipe-24.0-linux.zip

FRAGPIPE_BIN=/home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe
FRAGPIPE_TOOLS=/home/maropakis.a/fragpipe/fragpipe-24/tools

source ~/.bashrc
```

## FASTA File Generation   

FragPipe requires decoy sequences in the search database (target-decoy FDR estimation). You should maintain two separate copies of the FASTA:
* ```/scratch/maropakis.a/Dependencies/FASTA/``` -- FASTAs with no decoys OR ```/scratch/maropakis.a/Dependencies/FASTA_appended/```-- FASTAs with no decoys, custom sequences appended 
* ```/scratch/maropakis.a/Dependencies/FASTA_fragpipe/``` -- generated FASTAs with reversed decoy entries appended (```rev_``` prefixes). Decoys are full sequence reversals (not shuffled).

To generate the FragPipe FASTA copies, run ```1_PrepFASTA.py``` and ```2_genFragFASTA.py``` by submitting ```_build_fragFASTA.sh```. 
- 1_PrepFASTA.py -- creates species-specific `.csv` files used by `2_buildFragFASTA.py`. MTP peptides from appended FASTAs are matched to base peptides identified in MaxQuant searches.
- 2_buildFragFASTA.py -- Builds FragPipe-ready FASTAs from appended MTP FASTAs and the metadata tables generated in Step 1. Reference entries pass through unchanged. Kept MTP entries receive FragPipe-compatible headers and are added to the database. Reverse-sequence decoys are appended for every target entry.

>[!WARNING]
> This pipeline is optimized for the AAS Pipeline described in [Tsour et al., 2026](https://pubmed.ncbi.nlm.nih.gov/39253435/). If not intending to append custom peptides to the FASTA, skip ```1_PrepFASTA.py```. 

Inputs from MaxQuant Dependent Peptide Search: 
- evidence.txt
- proteinGroups.txt

Outputs:
- human.csv
- mouse.csv
- {species}.csv etc...

Example:
```bash
python3 1_PrepFASTA.py \
  --mtp-dir /scratch/maropakis.a/Dependencies/FASTA_appended/ \
  --human-root /scratch/maropakis.a/MQ_outputs/Ping_2018 \
  --mouse-root /scratch/maropakis.a/MQ_outputs/Takasugi_2024 \
  --out-dir /scratch/maropakis.a/Dependencies/mtp_maps/
```
```bash
python3 2_buildFragFASTA.py \
  --mtp-dir /scratch/maropakis.a/Dependencies/FASTA_appended/ \
  --human-csv /scratch/maropakis.a/Dependencies/mtp_maps/human.csv \
  --mouse-csv /scratch/maropakis.a/Dependencies/mtp_maps/mouse.csv \
  --out-dir /scratch/maropakis.a/Dependencies/FASTA_fragpipe/
```
or
```bash
sbatch _build_fragFASTA.sh
```

## Generating .workflow files 

The FragPipe ```.workflow``` file should be generated and exported from the FragPipe GUI. 

**Key dependencies may include:**
1. Search type: TMT closed search, standard 
2. TMT channels (TMT only)
3. TMT-integrator=(Philosopher for .mzML files, IonQuant for .raw files)
4. Decoy prefix: ```rev_``` (must match what ```gen_fragpipe_fasta.py``` uses)
5. Protein FDR: default template value (```--prot 0.05```)
6. Annotation file (for TMT searches only)

**To create the ```.workflow``` file:**
1. Open FragPipe GUI on a machine with test RAW files
2. Load RAW files --> select experiment type
3. Configure MSFragger search parameters for closed search
4. Enable Philosopher and TMT-integrator
5. Export workflow --> save as ```descriptivename.workflow```
6. Copy to cluster into ```templates``` folder

## Generating Annotation files 
```3_gen_annotations.py``` generates FragPipe annotation files from per-plex ```sample_map.xlsx``` files. 

Input: sample_map_*.xlsx
Output: *_annotation.txt 

## Staging spectra 
```4_stage_spectra.py``` builds one folder per-plex and stages:
- RAW files
- mzML files
- annotation.txt files
using symbolic links (symlinks).

Example:
```bash
python3 4_stage_spectra.py \
  --raw-root /scratch/maropakis.a/MQ_raw \
  --annot-dir /scratch/maropakis.a/Dependencies/annotations \
  --spectra-root /scratch/maropakis.a/spectra
```

Current dataset layouts handled by the script include:
```text
Ping_2018
  ├── ACG
  └── FC

Takasugi_2024
  ├── aorta
  ├── brain
  ├── heart
  ├── kidney
  ├── liver
  ├── lung
  ├── muscle
  └── skin
```

## Optional unless running Philosopher: Convert RAW files to mzML
```5_msconvert.sh``` uses ThermoRawFileParser to convert Thermo `.raw` files into `.mzML`. Submit as a SLURM array job.

Example:
```bash
sbatch --array=1-18 5_msconvert.sh
```

## Run FragPipe! 
```7_run_plexes.py``` creates workflow and manifest files for each plex and launches FragPipe headlessly. Consumes the output of `4_stage_spectra.py`.

Workflow routing:
| Dataset | Workflow |
|----------|----------|
| ACG / FC | TMT10_MS3_Val.workflow |
| aorta, brain, heart, kidney, liver, lung, muscle, skin | TMT16_Val.workflow |

Example:
```bash
python3 7_run_plexes.py \
  --spectra-root /scratch/maropakis.a/spectra \
  --fasta-dir /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
  --template-dir /home/maropakis.a/scripts/FragPipe/templates \
  --out-dir /scratch/maropakis.a/Frag_outputs \
  --fragpipe-bin /home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe \
  --tools-folder /home/maropakis.a/fragpipe/fragpipe-24.0/tools
```

Run a single plex:
```bash
python3 7_run_plexes.py --only acgb1 --run
```

or ```sbatch _submit_fragpipe.sh```. 


