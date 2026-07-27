# pipeline/

![Python](https://img.shields.io/badge/Python-3-3776AB?style=flat-square&logo=python&logoColor=white)
![FragPipe](https://img.shields.io/badge/FragPipe-24.0-6E44FF?style=flat-square)
![Scheduler](https://img.shields.io/badge/scheduler-SLURM-00838F?style=flat-square)
![Labelling](https://img.shields.io/badge/TMT-10%20%7C%2011%20%7C%2016-F9A825?style=flat-square)
![Modes](https://img.shields.io/badge/modes-per--plex%20%7C%20experiment--level-2E7D32?style=flat-square)

The scripts that turn raw TMT data into finished FragPipe searches. Three search modes plus
the shared setup tools they all depend on.

Everything here assumes the per-plex FragPipe FASTAs already exist — build those first with
`utility/` (see below).

## Layout

| Directory | What it does |
|---|---|
| `utility/` | Shared setup, run before anything else: build the per-plex FragPipe FASTAs (`prepFASTA.py` → `buildFragFASTA.py`) and, optionally, pre-convert raw → mzML (`msconvert.py`). |
| `per-plex/` | Search each TMT plex on its own — one workflow, one manifest, one submit script per plex. |
| `experiment-level-plex/` | Search several plexes together as one run with TMT bridge-channel normalization. |
| `label-free/` | Placeholder for a label-free variant of the pipeline (not yet implemented). |

## Typical order

1. **Build FASTAs** — `utility/build_fragFASTA.sh` (resolve SAAPs per plex, write search-ready
   FASTAs, check for duplicate headers).
2. **Prep searches** — run `per-plex/` or `experiment-level-plex/` depending on whether plexes
   are searched alone or bridged together. Each writes a workflow, manifest, annotation, and
   SLURM submit script into your `Frag_outputs` root.
3. **Submit** — `sbatch` the generated `submit_*.sh` scripts.

## Shared conventions

- **TMT only.** The per-plex and experiment-level generators are for TMT-labelled data; channel
  order is the canonical TMT10–TMTpro18 sequence.
- **Species** is `human` or `mouse`, used to tag headers (OS/OX) — not guessed from filenames.
- **Plex/run tokens** are lowercased but keep underscores (`cortex_keele`), and are the key that
  ties a sample_map, a FASTA, and a spectra directory together.
- **Outputs** land under one `--out-dir` (`Frag_outputs/`) split into `workflows/`, `manifests/`,
  `annotations/`, `submit/`, `results/`, and `logs/`.
- Paths in the `_run_.sh` examples point at one specific cluster account — edit them for your own.
