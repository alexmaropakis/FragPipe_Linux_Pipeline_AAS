# Per-Plex FragPipe FASTA Pipeline

A two-stage pipeline that builds **per-batch, search-engine-ready FASTAs** containing a set of
"target" sequences of interest (e.g. variant/modified peptides) plus their reference proteome,
resolved strictly against each batch's own upstream search evidence — then hands the result off
to a per-batch FragPipe (or similar) search.

This README describes the general pattern so it can be adapted to a different proteomics
project, target-sequence type, or search engine. The original implementation was built for
TMT-plex mistranslated-peptide (MTP/SAAP) detection, but the architecture generalizes to any
"set of candidate sequences that need per-batch provenance resolution before being added to a
search database."

## What problem this solves

If you have:
- multiple independent acquisition batches (plexes, runs, samples — whatever your unit of
  analysis is), each searched separately upstream, and
- a list of candidate/non-standard sequences (variants, modified peptides, novel ORFs, etc.)
  that need to be matched back to a parent/reference identity using **only that batch's own**
  upstream evidence,

...then pooling all batches together to resolve those candidates is a correctness bug: a
candidate seen in batch A should never be resolved using evidence from batch B, even if the
same sequence appears in both. This pipeline enforces that boundary structurally, not just by
convention.

## Pipeline shape

```
Stage 1 (resolve)        Stage 2 (build)              Stage 3 (search)
candidate seqs    --->   per-batch CSV     --->        per-batch FASTA   --->   per-batch search
+ per-batch            (sequence -> parent           (filtered + decoy-          job (your engine)
upstream evidence       accession/gene/desc)           appended, search-
                                                        engine-safe headers)
```

| Generic step | Original filenames | Purpose |
|---|---|---|
| Stage 1 — resolve | `1_PrepFASTA.py` | Per batch: match each candidate sequence to a parent identity using only that batch's own upstream evidence. Emit one CSV per batch. |
| Stage 2 — build | `2_buildFragFASTA.py` | Per batch: filter the candidate FASTA down to "resolved" entries using that batch's CSV, rewrite headers to be search-engine-safe, append decoys. |
| Driver | `build_FragFASTA.sh` | Runs Stage 1 → Stage 2 → a sanity check (e.g. no duplicate headers) as one SLURM/batch job. |
| Downstream orchestration | `run_all_plexes.sh` | The single place that enumerates every batch and kicks off its actual search job using the Stage-2 FASTA. |

## Core design principles (keep these when adapting)

1. **Strict per-batch resolution, no pooling.** Each batch's candidates are resolved only
   against that batch's own evidence directory. A batch whose token matches zero or more than
   one evidence directory should be a hard error, not a silent best-guess — misrouted evidence
   is worse than a crashed job.
2. **A single token scheme ties everything together.** One filename → token function (e.g.
   "strip known prefix/suffix, lowercase") is used to match: candidate-sequence files ↔
   evidence directories ↔ per-batch CSVs ↔ output FASTAs. Keep this logic in exactly one place
   per stage and make Stage 2 reuse Stage 1's token function (or an identical copy) — drift
   between the two is a common source of "silently dropped batch" bugs.
3. **Disambiguate collisions explicitly, don't guess.** If your batch/sample names repeat
   across source datasets or experiments (e.g. the same tissue name used in two different
   studies), keep an explicit lookup table (dataset → suffix) rather than relying on directory
   structure alone. Adding a new ambiguous dataset should be a one-line addition to that table.
4. **Species/condition/other metadata travels with the row, not the filename.** Anything Stage
   2 needs (species, tag, group) should be read from the per-batch CSV that Stage 1 wrote, not
   re-derived from the filename — that way Stage 2 needs zero edits when a new dataset is added.
5. **Decoys and header format match your search engine's parser exactly.** Search engines like
   Philosopher/MSFragger are strict and silently drop malformed entries rather than erroring —
   know your engine's required header format and accession-uniqueness rules before writing
   Stage 2, and add a uniqueness pass for accessions if a candidate ID could repeat.
6. **Loud failure over silent fallback.** Missing CSV, missing evidence directory, ambiguous
   token, duplicate headers — all of these should stop the pipeline with a clear error, since
   the cost of a silently wrong FASTA (wrong evidence used) is much higher than a failed job.

## Adapting this pipeline to your own project

Work through these in order — most of the actual code (the matching, CSV I/O, decoy logic) can
stay close to the original; what changes is the small set of project-specific definitions.

### 1. Define your "batch" and your "candidate sequence"
What's your independent resolution unit (plex, run, sample, patient)? What's the list of
non-standard sequences you need to resolve (variant peptides, isoforms, novel transcripts)?
This determines what Stage 1 reads as input.

### 2. Define your token function
Write one function that turns a candidate-file name into a canonical batch token, and a second
that turns an evidence-directory name into the same token space. These must agree by
construction — test them against a few real filenames before running anything else.

### 3. Define your disambiguation table
If two different source datasets can produce the same bare token (e.g. same tissue/sample name
used twice), add an explicit `DATASET_SUFFIX`-style dict keyed on dataset name. Don't try to
infer disambiguation from path depth or ordering — it's fragile.

### 4. Point Stage 1 at your upstream search evidence format
The original reads MaxQuant `evidence.txt` (sequence → protein string) and `proteinGroups.txt`
(accession → gene/description via `Fasta headers`). If your upstream engine is different
(e.g. FragPipe `psm.tsv`, Comet/Percolator output, a custom search), swap in the equivalent
sequence→parent-protein and accession→gene/description lookups, but keep the per-batch
indexing (`{token: evidence_dir}`) and per-batch column schema (`sequence, accession, gene,
description, status, ...`) the same shape so Stage 2 doesn't need to change.

### 5. Define your candidate→parent matching rule
The original matches by "same length, exactly one residue different" (a single amino-acid
substitution). Replace this with whatever matching rule fits your candidate type — e.g. exact
substring, edit distance ≤ k, shared peptide backbone with a PTM mass offset, etc. Keep it
isolated in one function so it's easy to swap.

### 6. Define your header rewrite and decoy convention
Write the mock header format your search engine expects, and decide whether decoys are
reversed, shuffled, or generated by your engine itself (in which case Stage 2 may not need to
add them at all — check whether your search engine's own decoy generation is being used
downstream, in which case appending decoys yourself causes double-decoying).

### 7. Write your driver and orchestration scripts last
Once Stages 1–2 work for one batch end-to-end, write the `build_*.sh` driver (resolve → build →
sanity check) and the `run_all_*.sh` enumeration script (one line per batch, pointing at the
finished FASTA + your search engine's workflow/config for that batch type). Keep the
"add a dataset" instructions co-located in a comment block at the top of the enumeration script
— it's the part that changes most often.

## Suggested sanity checks to keep regardless of project

- No duplicate headers in any output FASTA (a duplicate accession will cause silent collapses
  or ambiguous PSM assignment downstream).
- Every batch's candidate file matches **exactly one** evidence directory — zero or multiple
  matches should hard-fail.
- Spot-check that `status=keep` counts look reasonable per batch (a batch resolving 0% or 100%
  of candidates usually means a token mismatch, not a biological result).
