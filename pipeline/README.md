# Pipeline for Running FragPipe 24.0 in Linux

## Templates

| Template | Dataset | Species | Acquisition | TMT | Reporter level |
|---|---|---|---|---|---|
| `TMT10_MS3_Val.workflow` | Ping 2018 (ACG/FC) | human | MS3 | TMT-10 | 3 |
| `TMT16_Val.workflow` | Takasugi 2024 (tissues) | mouse | MS2 | TMT-16 | 2 |

Both are closed searches with MSBooster + Percolator, pointed at the
`_fragpipe.fasta` databases (reference + kept MTPs + `rev_` decoys) from the
BLASTp stage. Plexes route to a template by name token: `FC*`/`ACG*` → human/MS3;
tissue names (aorta, brain, heart, kidney, liver, lung, muscle, skin) → mouse/MS2.

## Pipeline order

```
build_fragpipe_fasta.py   (databases — see BLASTp stage)
        │
gen_annotations.py        sample_map/*.xlsx → annotations/<slug>_annotation.txt
        │
stage_spectra.py          MQ_raw/ → spectra/<plex>/ (symlinks + annotation.txt)
        │
run_plexes.py (dry run)   verify routing, FASTA match, channel counts
        │
submit_fragpipe.sh        SLURM array, one task per plex
```

## Scripts

- **`gen_annotations.py`** — one FragPipe `annotation.txt` per TMT set from the
  `sample_map_*.xlsx` files. Sorts channels into canonical order and
  disambiguates labels that repeat within a plex (e.g. GIS on two channels).
  Paths are set at the top of the file (`MAP_DIR`, `OUT_DIR`).
- **`stage_spectra.py`** — builds one folder per plex under `--spectra-root`,
  symlinking the spectra and copying in the matching `annotation.txt`. Handles
  flat Takasugi tissue folders and nested Ping ACG/FC batch subfolders (`ACG/b1`
  → `acgb1`, etc.).
- **`run_plexes.py`** — for each plex folder: routes to a template, patches
  `database.db-path` and `tmtintegrator.channel_num`, writes the `.workflow` and
  `.fp-manifest`, and (with `--run`) launches FragPipe headless. Without `--run`
  it only prints the commands — this is the dry run.
- **`submit_fragpipe.sh`** — SLURM array wrapper; each task reads one plex name
  from `plex_list.txt` and calls `run_plexes.py --only <plex> --run`.

## Usage

### 1. Generate annotations

```bash
python3 gen_annotations.py
```

Reads `sample_map/*.xlsx`, writes `annotations/<slug>_annotation.txt`. Each
`.xlsx` needs `tmt_channel` and `sample_name` columns (header casing/spacing is
normalized). The slug comes from the filename: `sample_map_acgb5.xlsx` → `acgb5`.

### 2. Stage spectra

```bash
python3 stage_spectra.py \
  --raw-root     /scratch/maropakis.a/MQ_raw \
  --annot-dir    /scratch/maropakis.a/Dependencies/annotations \
  --spectra-root /scratch/maropakis.a/spectra
```

The spectra root is *created* here — it's not a pre-existing directory. Each
resulting `spectra/<plex>/` holds symlinked spectra plus `annotation.txt`. A
`WARN ... NO annotation` line means the slug from `gen_annotations.py` didn't
match the plex name — fix the `sample_map` filename or the staging slug.

### 3. Dry run

```bash
bash run_dry.sh
```

which is:

```bash
python3 run_plexes.py \
  --spectra-root /scratch/maropakis.a/spectra \
  --fasta-dir    /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
  --template-dir /home/maropakis.a/scripts/FragPipe/templates \
  --out-dir      /scratch/maropakis.a/Frag_outputs \
  --fragpipe-bin /home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe \
  --tools-folder /home/maropakis.a/fragpipe/fragpipe-24.0/tools
```

No `--run`, so it only prints. Confirm every plex shows the right template,
channel count, and a matched `fasta=`. Investigate any `SKIP` (missing FASTA /
annotation / spectra, or no routing rule) or `WARN` (channel count ≠ expected)
before submitting. Add `--only <plex>` to check a single plex.

### 4. Submit the array

Build the plex list (one plex name per line, matching the staged folder names),
then submit:

```bash
ls /scratch/maropakis.a/spectra > /scratch/maropakis.a/Frag_outputs/plex_list.txt
sbatch --array=1-18 submit_fragpipe.sh
```

`--array=1-18` with no `%N` cap submits all 18 plexes as independent tasks with
no concurrency limit; add e.g. `%5` to cap at five concurrent. The array index
selects the line from `plex_list.txt`. Test one task first:

```bash
sbatch --array=1-1 submit_fragpipe.sh
```

`submit_fragpipe.sh` itself can be regenerated from the heredoc in
`submit_prep.txt` (copy/paste into the terminal).

## Editing the templates

`run_plexes.py` patches only `database.db-path` and `tmtintegrator.channel_num`
per plex. Everything else (decoy tag, MSBooster, extraction tool, MSstats) lives
in the templates and is set once. Edits we've made:

**Reporter extraction Philosopher → IonQuant** — Philosopher can't read `.raw`
(mzML only); IonQuant reads `.raw` natively. No `^...$` anchors (trailing
whitespace / CR makes anchored patterns miss silently):

```bash
sed -i.bak 's/extraction_tool=Philosopher/extraction_tool=IonQuant/' \
  templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
grep -H 'extraction_tool=' templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
```

**Disable MSstats** — `philosopher-msstats=true` requires Philosopher as the
extraction tool, so it conflicts with IonQuant. The output isn't used downstream
(`Validation2`/`Quant_TMT_SAAPs` consume `tmt-report` tables and `psm.tsv`):

```bash
sed -i.bak2 's/philosopher-msstats=true/philosopher-msstats=false/' \
  templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
grep -H 'philosopher-msstats=' templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
```

> `sed -i.bak` saves the original as `<file>.bak`; restore with
> `mv <file>.bak <file>`. The database already carries `rev_` decoys, so leave
> any add-decoys / add-contaminants option **off** (re-decoying breaks FDR).

### Edits don't propagate to already-generated per-plex copies

Editing a template does **not** update `.workflow` files `run_plexes.py` already
wrote into `Frag_outputs/`. Regenerate by re-running the dry run / submission.
If a run throws an error you thought you fixed, check the copy it actually used:

```bash
grep -r '^tmtintegrator.extraction_tool=' \
  /scratch/maropakis.a/Frag_outputs/results/*/ \
  /scratch/maropakis.a/Frag_outputs/workflows/
```

## Notes

- Java: `submit_fragpipe.sh` pins JDK 17 (`$HOME/bin/jdk-17.0.18+8`).
- The on-disk `submit_fragpipe.sh` requests 10 CPUs / 40G; the `submit_prep.txt`
  heredoc requests 16 CPUs / 64G. `run_plexes.py` defaults to `--threads 16
  --ram 64`, so the 10/40 version under-provisions threads relative to what
  FragPipe is told to use — pick one and keep them consistent.
- MSBooster needs a DIA-NN binary; confirm one exists or it errors mid-run:
  `find ~/fragpipe/fragpipe-24.0 -iname '*diann*' 2>/dev/null`.
- FragPipe output (`psm.tsv`/`ion.tsv`/`tmt-report`) is not drop-in compatible
  with the MaxQuant-based downstream scripts — a conversion step is still needed.
```
