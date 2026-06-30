# BrainDecode FASTA / Plex-Prep Pipeline

Builds per-plex, FragPipe-ready FASTAs (reference + MTP entries + reversed decoys) from
MaxQuant DP-search evidence, then preps each TMT plex for a FragPipe headless search.

Two stages, run in order:

1. **`1_PrepFASTA.py`** — resolve each MTP (mistranslated peptide) to a base-peptide accession,
   strictly within its own plex's MaxQuant DP evidence. Emits one CSV per plex.
2. **`2_buildFragFASTA.py`** — use those per-plex CSVs to filter the appended `*_MTP.fasta`
   files into FragPipe-safe FASTAs with Philosopher-compatible headers and decoys.

`build_FragFASTA.sh` runs both stages plus a duplicate-header sanity check. `run_all_plexes.sh`
is a separate, downstream script that takes the finished FragPipe FASTAs and generates the
per-plex FragPipe workflow/manifest/submit scripts for every TMT plex across all datasets.

## Core principle: strict per-plex resolution

Every plex (TMT batch / tissue) is resolved **only** against its own MaxQuant DP
`combined/txt` evidence. An MTP sequence seen in multiple plexes is never pooled — it's
resolved independently in each plex and written to that plex's own CSV. This was a prior bug
(cross-plex pooling produced incorrect on-disk FASTAs) and is now enforced end-to-end:
a filename matching zero or more than one DP dir is a hard error, not a silent fallback.

## Stage 1 — `1_PrepFASTA.py`

For each `*_MTP.fasta` in `--mtp-dir`:

- Derives a **plex token** from the filename (`S1_ACGB1_MTP.fasta` → `acgb1`).
- Matches that token to exactly one `*_DP/combined/txt/` directory found by walking the
  `--human-root` / `--mouse-root` paths (species comes from which root flag found the dir —
  adding a dataset means passing its root, no token rule to edit).
- For each MTP sequence, finds same-length `evidence.txt` peptides differing at exactly one
  residue (`find_base_peptides`), then resolves the first candidate with a real UniProt
  accession (`MTP|` entries excluded) via `proteinGroups.txt` → `(accession, gene, description)`.
- Writes one CSV per plex, `{token}.csv`, with columns:
  `sequence, species, accession, gene, description, bp_seq, all_accessions, status, n_base_candidates`
  (`status` = `keep` or `unresolved`).

**Token disambiguation:** tissue names repeat across studies (e.g. Takasugi kidney vs. Keele
kidney). Datasets listed in `DATASET_SUFFIX` (currently `keele2025`→`keele`,
`tsumagari_2023`→`tsumagari`) get a dataset suffix appended to their token
(`kidney`→`kidney_keele`, `cortex_1`→`cortex_1_tsumagari`). Datasets with globally-unique
plex names (Ping ACG/FC, Bai pooled, Takasugi tissues) keep the bare token. Adding a new
disambiguated dataset = one entry in `DATASET_SUFFIX`.

```
python3 1_PrepFASTA.py \
  --mtp-dir     /scratch/maropakis.a/Dependencies/FASTA_appended/ \
  --human-root  /scratch/maropakis.a/MQ_outputs/Ping_2018 \
  --human-root  /scratch/maropakis.a/MQ_outputs/Bai_2020 \
  --mouse-root  /scratch/maropakis.a/MQ_outputs/Takasugi_2024 \
  --mouse-root  /scratch/maropakis.a/MQ_outputs/Keele_2025 \
  --mouse-root  /scratch/maropakis.a/MQ_outputs/Tsumagari_2023 \
  --out-dir     /scratch/maropakis.a/Dependencies/mtp_maps/
```

## Stage 2 — `2_buildFragFASTA.py`

For each `*_MTP.fasta`, loads its matching `{token}.csv` (same token logic as Stage 1) and
keeps only MTP entries with `status == keep`:

- Reference entries pass through unchanged.
- Kept MTP entries get a Philosopher-safe mock-UniProt header:
  `>sp|{accession}-{MTPid}-{token}|{gene}-mut {gene} mistranslated {MTPid} OS=... OX=... GN={gene} PE=1 SV=1`.
  The accession is de-duplicated per file (`-d2`, `-d3`, ...) if it collides.
- `OS=`/`OX=` species tag is read from the CSV's `species` column — never guessed from the
  filename, so new datasets need no code change here.
- Reversed (`rev_`) decoys are appended for **every** target, including MTPs.

Header parsing note: the real MTP FASTA header is `>MTP|7998_0_base...` — the id regex is
`MTP\|(\d+)`, not `MTP\d+`, otherwise distinct MTPs collide once written to Philosopher.

```
python3 2_buildFragFASTA.py \
  --mtp-dir /scratch/maropakis.a/Dependencies/FASTA_appended/ \
  --csv-dir /scratch/maropakis.a/Dependencies/mtp_maps/ \
  --out-dir /scratch/maropakis.a/Dependencies/FASTA_fragpipe/
```

## Driver — `build_FragFASTA.sh`

SLURM `sbatch` script (`short` partition, 16 cpus, 32G, 2h). Runs Stage 1 → Stage 2 → a
duplicate-header sanity check over every `*_fragpipe.fasta` in the output dir, and exits
non-zero if any duplicates are found.

```
sbatch build_FragFASTA.sh
```

> Note: this script `cd`s to `/home/maropakis.a/scripts/search_gen/` before invoking
> `1_PrepFASTA.py` / `2_buildFragFASTA.py` by relative path — confirm both scripts actually
> live in `search_gen/` (vs. directly under `scripts/`) before submitting.

## Downstream — `run_all_plexes.sh`

Separate SLURM script; the single place listing every TMT plex across all datasets. For each
plex it calls `gen_fragpipe_plex.py` (`Alex_gen_fragpipe.py` at
`/home/maropakis.a/scripts/search_gen/`) to stage raw→mzML conversion, spectra, and write a
per-plex FragPipe workflow + manifest + `submit_<plex>.sh`. Channel count is auto-derived from
each sample map.

- **Live now:** Ping_2018 (10 plexes, TMT10 MS3) and Bai_2020 (pooled, TMT10 MS2), Takasugi_2024
  (8 tissues, TMT16 MS2).
- **Commented out / blocked on TODOs:** Keele_2025 and Tsumagari_2023 — need real acquisition
  workflow templates (`TODO_KEELE.workflow`, `TODO_TSUMAGARI.workflow`) and their FragPipe
  FASTAs built via Stages 1–2 first.
- **Excluded entirely** (different, DIA pipeline): `PD_2026`, `Giansanti_2022`.

```
sbatch run_all_plexes.sh
for s in /scratch/maropakis.a/Frag_outputs/submit/submit_*.sh; do sbatch "$s"; done
```

## Adding a new TMT dataset

1. Place its `*_MTP.fasta` in `FASTA_appended/`, named `S?_<token>_MTP.fasta`.
2. If its plex/tissue names collide with an existing dataset, add a `DATASET_SUFFIX` entry in
   `1_PrepFASTA.py`.
3. Pass its MaxQuant DP root via `--human-root`/`--mouse-root` and rerun `build_FragFASTA.sh`.
4. Add one `gen` line per plex to `run_all_plexes.sh`, pointing at the dataset's real
   acquisition workflow template and sample map.
