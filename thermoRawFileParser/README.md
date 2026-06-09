# RAW → mzML → FragPipe Pipeline

Converts Thermo `.raw` to mzML, then runs per-plex FragPipe headless searches
across all TMT plexes on the Explorer cluster (MGHPCC).

## Layout

```
$SPECTRA/<plex>/         # raw + mzML + annotation.txt, one folder per plex
$SPECTRA/logs/           # conversion logs (not a plex)
$DEP/FASTA_fragpipe/     # per-plex search FASTAs
$OUT/{logs,workflows,manifests,results}/
$SCRIPTS/{run_plexes.py,submit_fragpipe.sh,templates/}
```

Plexes route by name: tokens starting `ACG`/`FC` → human (TMT-10, MS3);
tissue names (aorta, brain, heart, …) → mouse (TMTpro-16, MS2).

## Prerequisites

- ThermoRawFileParser (self-contained Linux build) — no Mono required
- FragPipe 24.0 with bundled `tools/`; MSFragger jar stored separately
- Java 17 (for FragPipe); .NET only if using a framework-based TRFP build

---

## Step 1 — Convert RAW → mzML

SLURM array, one task per plex. Output is **plain** `.mzML` (`-f=2`, no `-g`):
gzipped mzML is not searchable by MSFragger, and skipping gzip avoids a separate
decompression step.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=raw2mzml
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=$SPECTRA/logs/raw2mzml_%A_%a.log
set -uo pipefail

TRFP=~/thermoRawFileParser/ThermoRawFileParser
SPECTRA=/path/to/spectra

# Exclude logs/ so array indices map cleanly to real plexes.
mapfile -t PLEXES < <(find "$SPECTRA" -mindepth 1 -maxdepth 1 -type d ! -name logs | sort)
PLEX="${PLEXES[$((SLURM_ARRAY_TASK_ID - 1))]}"

for raw in "$PLEX"/*.raw; do
    [ -e "$raw" ] || continue
    base=$(basename "$raw" .raw)
    out="$PLEX/$base.mzML"
    [ -s "$out" ] && { echo "SKIP $base"; continue; }
    "$TRFP" -i="$raw" -o="$PLEX" -f=2 -l=3
    [ -s "$out" ] && echo "OK $base" || echo "FAIL $base"
done
```

Set `--array=1-N` where N is the number of plex folders (exclude `logs/`), then:

```bash
sbatch msconvert.sh
```

### Verify before deleting any raw

```bash
find $SPECTRA -name '*.mzML.gz' | wc -l        # expect 0
find $SPECTRA -name '*.mzML'    | wc -l         # == raw count when complete
find $SPECTRA -name '*.raw'     | wc -l

# raw files with no matching mzML (re-convert these; should be empty):
comm -23 \
  <(find $SPECTRA -name '*.raw'  | sed 's/\.raw$//'  | sort) \
  <(find $SPECTRA -name '*.mzML' | sed 's/\.mzML$//' | sort)
```

Once mzML and raw counts match:

```bash
find $SPECTRA -name '*.raw' -delete
```

---

## Step 2 — Run FragPipe (per plex)

`run_plexes.py` builds a per-plex `.workflow` and `.fp-manifest`, then runs
FragPipe headless. Per plex it injects `database.db-path`,
`fragger.fragger-path`, and `tmtintegrator.channel_num`.

> MSFragger is **not** discovered via `--config-tools-folder` and putting its jar
> in `tools/` does nothing — its path must be set in each workflow via
> `fragger.fragger-path`. Pass it with `--msfragger-path`. IonQuant only needs the
> same treatment if a template sets `ionquant.run-ionquant=true`.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=fragpipe
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=$OUT/logs/fp_%A_%a.out
#SBATCH --error=$OUT/logs/fp_%A_%a.err
set -euo pipefail
export JAVA_HOME=$HOME/bin/jdk-17.0.18+8
export PATH=$JAVA_HOME/bin:$PATH

PLEX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $OUT/plex_list.txt)

python3 $SCRIPTS/run_plexes.py \
  --spectra-root   $SPECTRA \
  --fasta-dir      $DEP/FASTA_fragpipe \
  --template-dir   $SCRIPTS/templates \
  --out-dir        $OUT \
  --fragpipe-bin   /path/to/fragpipe-24.0/bin/fragpipe \
  --tools-folder   /path/to/fragpipe-24.0/tools \
  --msfragger-path /path/to/MSFragger-x.y.z.jar \
  --only "$PLEX" --run
```

Test one plex, then run the full array:

```bash
sbatch --array=1-1 submit_fragpipe.sh
sbatch --array=1-N submit_fragpipe.sh
```

`plex_list.txt`: one plex id per line, no `logs` entry.

---

## Notes

- Keep mzML uncompressed (`-f=2`, no `-g`) — MSFragger can't read `.mzML.gz`.
- Never delete `.raw` until mzML/raw counts match.
- `logs/` lives under `$SPECTRA` but is not a plex — exclude it from any
  directory-derived plex list, or array indexing will be off by one and the
  last plex will be skipped.
- Clear `results/<plex>/` (stale `.cal`/`.index`, prior `combined/`) before
  re-running a plex.
- mzML ≈ 1 GB/file — check `df -h /scratch` before bulk conversion.
```
