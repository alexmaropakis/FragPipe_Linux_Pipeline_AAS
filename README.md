# FragPipe_Linux_Tutorial: An updated guide for running FragPipe in Linux.

[FragPipe](https://fragpipe.nesvilab.org/) is a comprehensive computational platform designed to analyze mass spectrometry-based proteomics data. It is powered by [MSFragger](https://msfragger.nesvilab.org/), an ultrafast proteomic search engine for both closed and open (wide precursor mass tolerance) peptide identification. It includes Percolator and the [Philosopher](https://nesvilab.github.io/philosopher/) toolkit for downstream post-processing of MSFragger search results, FDR filtering, and multi-experiment summary report generation. It also includes the [MSBooster](https://www.nature.com/articles/s41467-023-40129-9) module for deep-learning based rescoring of peptide identifications, and [Crystal-C](https://www.nesvilab.org/Crystal-C/) and [PTM-shepherd](https://github.com/Nesvilab/PTM-Shepherd) to aid interpretation of results from "open" and "mass offset" searches for PTMs.

This guide is targeted toward beginners planning to run FragPipe v24.0. It is optimized for [Northeastern University's Research Computing Cluster](https://rc.northeastern.edu/), which is hosted at the [Massachusetts Green High-Performance Computing Center (MGHPCC)](https://mghpcc.org/) in Holyoke, MA. 

Last updated 06-03-2026.

## Pipeline
1. Download FragPipe version into home directory
2. Download software framework dependencies
3. Create template workflow files in the FragPipe GUI
4. Move FASTA files into Linux instance and generate FASTAs with decoys (and MTPs appended if needed, as per Tsour et al. (2026) at the Slavov Lab)
5. Generate workflow file with ```gen_workflow.py```.
6. Run FragPipe!

## Setting up FragPipe in Linux

The most recent FragPipe version (and older versions) can be accessed via an academic license [here](https://github.com/Nesvilab/FragPipe/releases). The official tutorial for using FragPipe can be found [here](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html).

Download the wanted zip file and put it in your home directory.

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
wget https://openjdk-sources.osci.io/openjdk17/openjdk-17.0.18+8.tar.xz
tar -xvf openjdk-17.0.18+8.tar.xz
./configure
make
sudo make install
JAVA_HOME=/home/maropakis.a/bin/jdk-17.0.18+8
source ~/.bashrc
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

## FASTA Database Structure 

FragPipe requires decoy sequences in the search database (target-decoy FDR estimation). You should maintain two separate copies of the FASTA:
* ```/scratch/maropakis.a/Dependencies/FASTA/``` -- FASTAs with no decoys
* ```/scratch/maropakis.a/Dependencies/FASTA_fragpipe/``` -- FASTAs with reversed decoy entries appended (```rev_``` prefixes). Decoys are full sequence reversals (not shuffled).

To generate the FragPipe FASTA copies, run ```generation_scripts/gen_fragpipe_fasta.py```. 

## Generating .workflow files 

The FragPipe ```.workflow``` file should be generated and exported from the FragPipe GUI. 

**Key dependencies may include: **
1. Search type: TMT closed search, standard 
2. TMT channels (TMT only)
3. Decoy prefix: ```rev_``` (must match what ```gen_fragpipe_fasta.py``` uses)
4. Protein FDR: default template value (```--prot 0.05```)
5. Annotation file (for TMT searches only)

**To create the ```.workflow``` file: **
1. Open FragPipe GUI on a machine with test RAW files
2. Load RAW files --> select experiment type
3. Configure MSFragger search parameters for closed search
4. Enable Philosopher and TMT-integrator
5. Export workflow --> save as ```descriptivename.workflow```
6. Copy to cluster and pass to ```run_ms_search.sh```

## Run search 
1. Set up manifests + submit jobs by running ```sbatch run_ms_search.sh```
2. OR generate manifests first, add workflow manually, then submit




## 

