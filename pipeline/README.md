# pipeline/

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
