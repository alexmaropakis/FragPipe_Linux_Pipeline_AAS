# per-plex/

Prep one TMT plex for FragPipe, end to end, and get a submit script back. Use this when each
plex is searched and quantified on its own (no cross-plex bridging).

## Files

- `gen_fragpipe_plex.py` — the generator. TMT-labelled data only.
- `_run_.sh` — a SLURM wrapper holding a `gen` helper; list one `gen` call per plex here instead
  of typing them at the terminal.

## What the generator does, per plex

1. **annotation.txt** from the plex's `sample_map.xlsx` — one `<channel> <sample>` line, ordered
   by TMT channel. A label repeated within the plex (e.g. a reference on two channels) is made
   unique with a `_<channel>` suffix.
2. **raw → mzML** — converts every `.raw` in the raw dir with ThermoRawFileParser (plain indexed
   mzML), skipping files already converted. Skip entirely with `--no-convert`.
3. **spectra staging** — builds `<spectra-root>/<plex>/` with symlinked mzML plus the annotation.
4. **workflow + manifest** — patches the chosen `.workflow` template (`database.db-path`,
   `tmtintegrator.channel_num`) and writes the `.fp-manifest`.
5. **submit script** — `submit_<plex>.sh`, a headless FragPipe SLURM job for this plex.

The FASTA is found automatically in `--fasta-dir` by matching the plex token against
`*_fragpipe.fasta` names (the FASTAs built by `utility/buildFragFASTA.py`).

## Usage

```bash
python3 gen_fragpipe_plex.py <raw_dir> \
  --plex        acgb1 \
  --species     human \
  --workflow    templates/TMT10_MS3_Val.workflow \
  --sample-map  sample_map/acgb1.xlsx \
  --fasta-dir   Dependencies/FASTA_fragpipe \
  --out-dir     Frag_outputs \
  --spectra-root spectra

sbatch Frag_outputs/submit/submit_acgb1.sh
```

`--channels` is optional — the count comes from the sample_map; pass it only to assert an
expected value. Submit-script knobs (`--fragpipe-bin`, `--tools-folder`, `--java-home`,
`--partition`, `--threads`, `--ram`, `--time`) all have defaults you'll likely need to edit for
your cluster.

## Batching many plexes

Add a `gen` line per plex in `_run_.sh`, then submit them together:

```bash
for s in "$OUT"/submit/submit_*.sh; do sbatch "$s"; done
```
