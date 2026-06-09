# Pipeline to run FragPipe 24.0 in Linux. 
Note: this tutorial is optimized in the context of research I am doing in the Slavov Lab studying substituted amino acid peptides (SAAPs) arising from alternative RNA decoding. Thus, there will be discussion of appending SAAP sequences to canonical .FASTA files with FragPipe-compatible headers. If not performing a similar analysis, run the Pipeline found in `Generic`. 

## Requirements

- NCBI BLAST+ 2.17.0 — `~/bin/ncbi-blast-2.17.0+/bin/` (must be on `PATH`)
- FragPipe 24.0 — `~/fragpipe/fragpipe-24.0/` (MSFragger, IonQuant, Philosopher,
  TMT-Integrator, MSBooster, Percolator, DIA-NN bundled)
- JDK 17 — `~/bin/jdk-17.0.18+8`
- Python 3 with `pandas`
- SLURM (`short` partition)

Batch jobs don't read `.bashrc` aliases, so the submit scripts export `PATH` /
`JAVA_HOME` explicitly.

---

## Stage A — BLASTp and database construction

### A1. Map SAAPs to parent proteins — `blast_mtps.py`

Reads `SAAP_quant_df.csv` (columns `MTP_seq`, `BP_seq`, `Dataset`), dedupes on
`(dataset, MTP, BP)`, computes the substitution (recorded **BP→MTP**), splits by
species, and BLASTs each unique SAAP against the matching reference proteome with
`blastp-short` (`-evalue 1000`, `-comp_based_stats 0`, `-max_target_seqs 50`).

A hit counts as full-length only with no gaps, query start at 1, and full query
coverage. Classification of each SAAP:

| status | meaning |
|---|---|
| `drop_reference` | exact 0-mismatch full-length match → it's just the reference |
| `keep` | best full-length 1-mismatch parent → genuine substitution |
| `drop_no_parent` | no 1-mismatch full-length parent |
| `drop_trypsin` | best parent is trypsin (`PRSS1/2/3`, `TRY*`) |
| `drop_ig` | best parent is immunoglobulin (`IGH/K/L`, `JCHAIN`, `IGJ`) |

`keep` rows are flagged `ambiguous` when 1-mismatch candidates span >1 gene. The
subject start is captured to report the absolute protein position of the
substitution (`protein_pos`).

Run (matches `blast_mtps.sh`):

```bash
sbatch blast_mtps.sh
```

which calls:

```bash
python3 blast_mtps.py \
  --csv        $DEP/mtp_maps/SAAP_quant_df.csv \
  --human-ref  $DEP/FASTA/HUMAN.fasta \
  --mouse-ref  $DEP/FASTA/MOUSE_UP000000589_10090.fasta \
  --out-dir    $DEP/mtp_maps/ \
  --prefix     April26 \
  --threads    16
```

Outputs per species: `<prefix>_{human,mouse}_unfiltered.csv` (all SAAPs with
classification), `<prefix>_{human,mouse}_filtered.csv` (`keep` only), and the
intermediate query FASTA + raw BLAST TSV. The reference DB is built on demand
with `makeblastdb -parse_seqids`.

> Column flags default to `MTP_seq`, `BP_seq`, and `\ufeffDataset` (the dataset
> column carries a BOM). Override with `--seq-col` / `--bp-col` / `--dataset-col`
> if the CSV headers differ; the script lists available columns on mismatch.

### A2. Build search databases — `build_fragpipe_fasta.py`

For each `*_MTP.fasta`: reference entries pass through unchanged; kept MTPs (matched
by sequence against the filtered CSV) get a Philosopher-safe mock-UniProt header
with a unique accession `<acc>-<mtp_id>-<tmt_set>`; reversed `rev_` decoys are
appended for **every** target. Output is `<name>_fragpipe.fasta`.

```bash
python3 build_fragpipe_fasta.py \
  --mtp-dir   $DEP/FASTA_appended/ \
  --human-csv $DEP/mtp_maps/human_filtered.csv \
  --mouse-csv $DEP/mtp_maps/mouse_filtered.csv \
  --out-dir   $DEP/FASTA_fragpipe/
```

`build_fragfasta.sh` chains A1 + A2 plus a duplicate-header sanity check in one
SLURM job. Because the FASTA already carries `rev_` decoys, FragPipe must **not**
add decoys or contaminants (re-decoying breaks FDR).

---

## Stage B — FragPipe search

### Templates

| Template | Dataset | Species | Acquisition | TMT | Reporter level |
|---|---|---|---|---|---|
| `TMT10_MS3_Val.workflow` | Ping 2018 (ACG/FC) | human | MS3 | TMT-10 | 3 |
| `TMT16_Val.workflow` | Takasugi 2024 (tissues) | mouse | MS2 | TMT-16 | 2 |

Both: closed search, MSBooster + Percolator on, IonQuant as the TMT extraction
tool, MSstats off, pointed at the `_fragpipe.fasta` databases. Plexes route by
name token (`FC*`/`ACG*` → human/MS3; tissue names → mouse/MS2).

### B1. Generate annotations — `gen_annotations.py`

One FragPipe `annotation.txt` per TMT set from `sample_map_*.xlsx` (needs
`tmt_channel` + `sample_name` columns; casing/spacing normalized). Channels are
sorted into canonical order and labels repeating within a plex are
disambiguated (e.g. GIS on two channels). Slug comes from the filename
(`sample_map_acgb5.xlsx` → `acgb5`). Paths set at top of file.

```bash
python3 gen_annotations.py
```

### B2. Stage spectra — `stage_spectra.py`

Builds one folder per plex under `--spectra-root` (created here, not
pre-existing), symlinking spectra and copying in the matching `annotation.txt`.
Handles flat Takasugi tissue folders and nested Ping ACG/FC batch subfolders
(`ACG/b1` → `acgb1`).

```bash
python3 stage_spectra.py \
  --raw-root     /scratch/maropakis.a/MQ_raw \
  --annot-dir    /scratch/maropakis.a/Dependencies/annotations \
  --spectra-root /scratch/maropakis.a/spectra
```

A `WARN ... NO annotation` line means the annotation slug didn't match the plex
name — fix the `sample_map` filename or staging slug.

### B3. Dry run — `run_plexes.py`

Routes each plex folder to a template, patches `database.db-path` and
`tmtintegrator.channel_num`, writes the `.workflow` + `.fp-manifest`, and (with
`--run`) launches FragPipe headless. Without `--run` it only prints — the dry run:

```bash
bash run_dry.sh
```

Confirm every plex shows the right template, channel count, and a matched
`fasta=`. Investigate any `SKIP` (missing FASTA / annotation / spectra, or no
routing rule) or `WARN` (channel count ≠ expected) before submitting. Add
`--only <plex>` for a single plex.

### B4. Submit the array — `submit_fragpipe.sh`

Build the plex list (one name per line, matching staged folder names), then
submit one task per plex:

```bash
ls /scratch/maropakis.a/spectra > /scratch/maropakis.a/Frag_outputs/plex_list.txt
sbatch --array=1-18 submit_fragpipe.sh
```

`--array=1-18` with no `%N` cap submits all 18 with no concurrency limit (add
e.g. `%5` to cap). Test one task first with `--array=1-1`. The script itself can
be regenerated from the heredoc in `submit_prep.txt`.

---

## Editing the templates

`run_plexes.py` patches only `database.db-path` and `tmtintegrator.channel_num`
per plex; everything else lives in the templates and is set once. Edits made:

**Reporter extraction Philosopher → IonQuant** — Philosopher reads mzML only;
IonQuant reads `.raw` natively. No `^...$` anchors (trailing whitespace / CR makes
anchored patterns miss silently):

```bash
sed -i.bak 's/extraction_tool=Philosopher/extraction_tool=IonQuant/' \
  templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
grep -H 'extraction_tool=' templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
```

**Disable MSstats** — `philosopher-msstats=true` requires Philosopher as the
extraction tool, so it conflicts with IonQuant; the output isn't used downstream:

```bash
sed -i.bak2 's/philosopher-msstats=true/philosopher-msstats=false/' \
  templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
grep -H 'philosopher-msstats=' templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
```

> `sed -i.bak` saves the original as `<file>.bak`; restore with
> `mv <file>.bak <file>`.

**Edits don't propagate to already-generated per-plex copies.** Editing a
template does not update the `.workflow` files `run_plexes.py` already wrote into
`Frag_outputs/`; regenerate by re-running the dry run / submission. If a run
throws an error you thought you fixed, check the copy it actually used:

```bash
grep -r '^tmtintegrator.extraction_tool=' \
  /scratch/maropakis.a/Frag_outputs/results/*/ \
  /scratch/maropakis.a/Frag_outputs/workflows/
```

---

## Key paths

| | |
|---|---|
| Scripts | `/home/maropakis.a/scripts/` (`BLASTp/`, `FragPipe/`) |
| Dependencies | `/scratch/maropakis.a/Dependencies/` |
| Reference FASTAs | `Dependencies/FASTA/{HUMAN,MOUSE_UP000000589_10090}.fasta` |
| Appended MTP FASTAs | `Dependencies/FASTA_appended/*_MTP.fasta` |
| Search databases | `Dependencies/FASTA_fragpipe/*_fragpipe.fasta` |
| Annotations | `Dependencies/annotations/` |
| Staged spectra | `/scratch/maropakis.a/spectra/<plex>/` |
| FragPipe outputs | `/scratch/maropakis.a/Frag_outputs/` |
