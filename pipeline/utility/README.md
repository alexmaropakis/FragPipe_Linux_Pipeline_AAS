# utility/

![Stage](https://img.shields.io/badge/stage-run%20first-D84315?style=flat-square)
![Python](https://img.shields.io/badge/Python-3-3776AB?style=flat-square&logo=python&logoColor=white)
![Deps](https://img.shields.io/badge/deps-pandas%20%7C%20ThermoRawFileParser-795548?style=flat-square)
![Resolution](https://img.shields.io/badge/SAAP%20resolution-per--plex-2E7D32?style=flat-square)
![Decoys](https://img.shields.io/badge/decoys-rev__-455A64?style=flat-square)

Shared setup that the search generators depend on. Two independent jobs:

1. **Build the per-plex FragPipe FASTAs** (`prepFASTA.py` → `buildFragFASTA.py`). Run this before
   any search — the `per-plex/` and `experiment-level-plex/` generators both expect the
   `*_fragpipe.fasta` files these produce.
2. **Convert raw → mzML** (`msconvert.py`), optional/standalone — the search generators already
   convert on the fly, so this is only for pre-staging spectra separately.

## FASTA build (two stages + a check)

The FASTA build resolves candidate SAAP/MTP sequences to a parent protein **per plex**, using only
that plex's own MaxQuant dependent-peptide evidence — never pooled across plexes.

### `prepFASTA.py` — Stage 1 (resolve)

Per plex, matches each substituted sequence in `*_MTP.fasta` to a parent accession/gene/description
using that plex's MaxQuant output (`evidence.txt` for proteins, `proteinGroups.txt` for
accession → gene/description), and writes one CSV per plex.

- **Inputs:** `*_MTP.fasta` (reference + appended substituted sequences), plus MaxQuant
  dependent-peptide search outputs grouped under `--human-root` / `--mouse-root` (repeatable).
- **Output:** `{token}.csv` with columns `sequence, accession, gene, description, bp_seq,
  all_accessions, status, n_base_candidates`.
- **Prerequisite:** the SAAP_Detection & Validation steps from the Decode pipeline
  (Tsour et al., *Nature* 2026).

Plex tokens are derived from the FASTA filename (`S1_ACGB1_MTP.fasta` → `acgb1`). When the same
tissue name appears across studies, `DATASET_SUFFIX` disambiguates it — add a one-line entry there
for each new colliding dataset.

### `buildFragFASTA.py` — Stage 2 (build)

Per plex, filters `*_MTP.fasta` down to the sequences its CSV marks `status=keep`, rewrites headers
into the search-engine-safe form, and appends `rev_` decoys. Species (for the OS/OX header tag) is
read from the CSV row, not guessed from the filename.

- **Output:** `*_fragpipe.fasta` (same `*` stem as the matching `*_MTP.fasta`).

### `build_fragFASTA.sh` — driver

Runs Stage 1 → Stage 2 → a duplicate-header check as one SLURM job. The duplicate check hard-fails
the job if any output FASTA has repeated `>` headers (duplicates cause silent collapses downstream).
Edit the account paths and the `--*-root` dataset lists at the top before running.

```bash
sbatch build_fragFASTA.sh
```

## raw → mzML

### `msconvert.py`

Converts every `.raw` in a directory to plain indexed mzML with ThermoRawFileParser (skipping any
already converted), then symlinks the results into `<spectra-root>/<plex>/`. Plain mzML matters:
gzipped `.mzML.gz` is silently skipped by FragPipe.

```bash
python3 msconvert.py <raw_dir> --plex <token> --spectra-root spectra
```

Requires ThermoRawFileParser on `PATH` (or pass `--trfp`).

### `msconvert.sh`

A one-off SLURM wrapper around `msconvert.py` for a single plex — edit the path and plex, then
`sbatch`.
