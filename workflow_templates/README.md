# Workflows for FragPipe

FragPipe 24.0 headless workflow templates for the BrainDecode MTP/SAAP search
pipeline, plus the commands used to generate and edit them.

Two templates, one per acquisition type:

| Template | Dataset | Acquisition | TMT | Reporter quant level |
|---|---|---|---|---|
| `TMT10_MS3_Val.workflow` | Ping 2018 (human, ACG/FC) | MS3 | TMT10 | 3 |
| `TMT16_Val.workflow` | Takasugi 2024 (mouse tissues) | MS2 | TMTpro 16 | 2 |

Both are closed searches with MSBooster + Percolator enabled, pointed at the
`_fragpipe.fasta` databases (reference + MTP entries + `rev_` decoys) produced by
the BLASTp stage.

## Where these came from

The templates are edited copies of FragPipe's bundled workflows. List the
bundled set:

```bash
ls ~/fragpipe/fragpipe-24.0/workflows/
```

`TMT10_MS3_Val.workflow` started from `TMT10-MS3.workflow` (MS3 reporter
extraction + Percolator). `TMT16_Val.workflow` started from the plain
`TMT16.workflow` (MS2 reporter extraction). Choosing the right base template is
what fixes the reporter-extraction MS level — that key is not edited by hand.

## Required per-template edits

The base templates need these settings for the MTP pipeline. They're applied to
the templates once; `run_plexes.py` then patches the per-plex specifics
(`database.db-path`, `tmtintegrator.channel_num`) at submission time.

Database settings — the FASTA already carries `rev_` decoys, so do **not** enable
any add-decoys/add-contaminants option (re-decoying doubles them and breaks FDR):

```
database.decoy-tag=rev_
msbooster.run-msbooster=true
```

### TMT intensity extraction: Philosopher → IonQuant

Philosopher cannot read `.raw` files (mzML only); IonQuant reads `.raw` natively.
Switch both templates. Note: no `^...$` anchors — the lines can carry trailing
whitespace / carriage returns that make anchored patterns silently miss.

```bash
sed -i.bak 's/extraction_tool=Philosopher/extraction_tool=IonQuant/' \
  templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow

grep -H 'extraction_tool=' templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
```

Both should now read `IonQuant`.

### Disable MSstats generation

`tmtintegrator.philosopher-msstats=true` requires Philosopher as the extraction
tool, so it conflicts once you've moved to IonQuant. The MSstats output isn't
consumed downstream (`Validation2` / `Quant_TMT_SAAPs` use the `tmt-report`
tables and `psm.tsv`), so turn it off:

```bash
sed -i.bak2 's/philosopher-msstats=true/philosopher-msstats=false/' \
  templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow

grep -H 'philosopher-msstats=' templates/TMT16_Val.workflow templates/TMT10_MS3_Val.workflow
```

Both should read `false`.

> `sed -i.bak` writes the edited file and saves the untouched original as
> `<file>.bak`. Restore with `mv <file>.bak <file>` if an edit breaks something.
> Drop the suffix (`sed -i`) for no backup.

## Verifying the edits propagated

`run_plexes.py` writes a fresh per-plex workflow into each output dir at
submission time, generated from the template **as it was then**. Editing a
template does *not* retroactively touch already-generated per-plex copies — those
must be regenerated. If a run throws an error you thought you'd fixed, check the
copy the run actually used, not just the template:

```bash
grep -r '^tmtintegrator.extraction_tool=' \
  /scratch/maropakis.a/Frag_outputs/results/*/ \
  /scratch/maropakis.a/Frag_outputs/workflows/
```

If those still show the old value, re-run `run_plexes.py` to regenerate them.

## Running headless

`run_plexes.py` routes each plex to the correct template by name token (ACG/FC →
human/MS3; tissue names → mouse/MS2) and launches FragPipe. The underlying
headless invocation is:

```bash
~/fragpipe/fragpipe-24.0/bin/fragpipe --headless \
  --workflow  <template>.workflow \
  --manifest  <plex>.fp-manifest \
  --workdir   <output_dir> \
  --config-tools-folder ~/fragpipe/fragpipe-24.0/tools/ \
  --ram 0 --threads 16
```

MSBooster requires a DIA-NN binary for its MS/MS + RT prediction; confirm one
exists or MSBooster errors mid-run:

```bash
find ~/fragpipe/fragpipe-24.0 -iname '*diann*' -o -iname '*DIA-NN*' 2>/dev/null
```

## Files

- `TMT10_MS3_Val.workflow` — human / Ping 2018, MS3, TMT10
- `TMT16_Val.workflow` — mouse / Takasugi 2024, MS2, TMTpro 16
```
