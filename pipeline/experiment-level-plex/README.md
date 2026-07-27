# experiment-level-plex/

Prep a FragPipe run that spans several TMT plexes searched and quantified together, so
TMT-Integrator can bridge across plexes on a shared reference channel. Use this instead of
`per-plex/` whenever plexes need to be normalized against each other.

## Files

- `gen_fragpipe_experiment_plex.py` — the generator. TMT-labelled data only.
- `_run_.sh` — a SLURM wrapper; list the full commands here (a loop over tissues works well when
  naming is consistent).

## Vocabulary

- **run** — the whole FragPipe job: one workflow, one manifest, one FASTA, one submit script.
- **experiment** — one TMT plex inside that run. You name it (`cortex_1`, `cortex_2`, …), point it
  at that plex's raw dir, and give the plex key that selects its rows from the shared
  sample_map.

## Sample map layout

One `.xlsx` for the whole run, with a plex column (`TMT plex` or `sample_ID`) so each experiment's
rows can be pulled out by key. The key is matched case-insensitively, so both `3` and `S3` work.

## What the generator does

For each `--experiment`: builds its annotation.txt from its plex rows, converts/stages its spectra
into `<spectra-root>/<run>/<experiment>/`. Then, once per run: writes one manifest listing every
file across all experiments, patches the workflow template, writes a human-readable sample table,
and emits `submit_<run>.sh`.

All experiments must agree on channel count (one workflow → one `channel_num`); a mismatch warns
and uses the max.

## Usage

```bash
python3 gen_fragpipe_experiment_plex.py \
  --run        cortex_tsumagari \
  --species    mouse \
  --workflow   templates/TMT10_MS2_Val.workflow \
  --sample-map sample_map/sample_map_tsumagari_cortex.xlsx \
  --experiment cortex_1 MQ_raw/Tsumagari_2023/cortex_1 3 \
  --experiment cortex_2 MQ_raw/Tsumagari_2023/cortex_2 4 \
  --fasta      Dependencies/FASTA_fragpipe/S9_cortex_tsumagari_fragpipe.fasta \
  --out-dir    Frag_outputs \
  --spectra-root spectra

sbatch Frag_outputs/submit/submit_cortex_tsumagari.sh
```

Repeat `--experiment NAME RAW_DIR PLEX_KEY` once per plex in the run. Unlike `per-plex/`, the FASTA
is passed explicitly with `--fasta` rather than looked up by token. The same submit-script and
conversion knobs (`--no-convert`, `--threads`, `--ram`, etc.) apply.
