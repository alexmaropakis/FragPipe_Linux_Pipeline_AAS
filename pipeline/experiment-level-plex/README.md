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

Bridge channels get a unique annotation name across all plexes, sharing the workflow's
`tmtintegrator.ref_tag` as a prefix: any `sample_name` containing it (e.g. `Bridge`) becomes
`Bridge1`, `Bridge2`, ..., numbered in `--experiment` order and then channel order.

All experiments must agree on channel count (one workflow → one `channel_num`); a mismatch warns
and uses the max.

## Usage

### Required arguments only

```bash
args=(
  --run          cortex_tsumagari                                                   # names all outputs
  --workflow     /home/$USER/scripts/FragPipe/templates/TMT10_MS2_Val.workflow
  --sample-map   /scratch/$USER/Dependencies/sample_map/sample_map_tsumagari_cortex.xlsx
  --experiment   cortex_1 /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_1 3           # NAME | DIR | PLEX_KEY
  --experiment   cortex_2 /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_2 4           # one per plex
  --fasta        /scratch/$USER/Dependencies/FASTA_fragpipe/S9_cortex_tsumagari_fragpipe.fasta
  --out-dir      /scratch/$USER/Frag_outputs
  --spectra-root /scratch/$USER/spectra
)
python3 gen_fragpipe_experiment_plex.py "${args[@]}"
sbatch /scratch/$USER/Frag_outputs/submit/submit_cortex_tsumagari.sh
```

### Every argument

```bash
args=(
  --run          cortex_tsumagari                                    
  --workflow     /home/$USER/scripts/FragPipe/templates/TMT10_MS2_Val.workflow
  --sample-map   /scratch/$USER/Dependencies/sample_map/sample_map_tsumagari_cortex.xlsx
  --experiment   cortex_1 /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_1 3         
  --experiment   cortex_2 /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_2 4        
  --fasta        /scratch/$USER/Dependencies/FASTA_fragpipe/S9_cortex_tsumagari_fragpipe.fasta
  --out-dir      /scratch/$USER/Frag_outputs
  --spectra-root /scratch/$USER/spectra
  --trfp         ~/thermoRawFileParser/ThermoRawFileParser                          # .raw -> .mzML
  --msconvert    msconvert                                                          # .mzXML -> .mzML
  --no-convert                                                                      # skip .raw conversion
  --allow-raw                                                                       # keep .raw in manifest
  --fragpipe-bin /home/$USER/fragpipe/fragpipe-24.0/bin/fragpipe
  --tools-folder /home/$USER/fragpipe/fragpipe-24.0/tools
  --java-home    ~/bin/jdk-17.0.18+8
  --partition    short                                                
  --threads      16
  --ram          64                                                     
  --time         24:00:00                                                  
)
python3 gen_fragpipe_experiment_plex.py "${args[@]}"
```

Optional values shown are the defaults; `--no-convert` and `--allow-raw` are off unless passed.
