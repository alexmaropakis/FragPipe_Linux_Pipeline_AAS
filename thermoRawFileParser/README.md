# RAW → mzML → FragPipe Pipeline

Converts Thermo `.raw` to plain mzML, then runs per-plex FragPipe headless
searches across 18 TMT plexes on Explorer (MGHPCC).

## Layout

```
/scratch/maropakis.a/spectra/<plex>/        # .raw + .mzML + annotation.txt
/scratch/maropakis.a/spectra/logs/          # conversion logs (NOT a plex)
/scratch/maropakis.a/Dependencies/FASTA_fragpipe/   # per-plex FASTAs
/scratch/maropakis.a/Frag_outputs/{logs,workflows,manifests,results}/
/home/maropakis.a/scripts/FragPipe/{run_plexes.py,submit_fragpipe.sh,templates/}
~/thermoRawFileParser/ThermoRawFileParser   # self-contained Linux build
/home/maropakis.a/fragpipe/fragpipe-24.0/   # FragPipe + bundled tools/
/home/maropakis.a/fragpipe/MSFragger-4.4.1/MSFragger-4.4.1.jar
```

18 plexes: `acgb1–5`, `fcb1–5` (human, TMT-10/MS3) and
`aorta, brain, heart, kidney, liver, lung, muscle, skin` (mouse, TMTpro-16/MS2).

---

## Step 1 — Convert RAW → mzML

`msconvert.sh` (SLURM array, one task per plex). Writes **plain** `.mzML`
(`-f=2`, no `-g`) so files are searchable directly with no gunzip step.

> **Two bugs that bit earlier, now fixed below:**
> - `logs/` is a subdir of `spectra/` and was being counted as a 19th "plex",
>   shifting array indices and dropping `skin`. The plex list excludes `logs`.
> - The old script used `-g` (gzip) but checked for `$base.mzML`, so the
>   exists-check never matched and every file logged FAIL. Output is now plain
>   `.mzML` and the check matches.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=raw2mzml
#SBATCH --partition=short
#SBATCH --array=1-18
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/maropakis.a/spectra/logs/raw2mzml_%A_%a.log
set -uo pipefail

TRFP=~/thermoRawFileParser/ThermoRawFileParser
SPECTRA=/scratch/maropakis.a/spectra

# Exclude logs/ so indices align and skin (last plex) is scheduled.
mapfile -t PLEXES < <(find "$SPECTRA" -mindepth 1 -maxdepth 1 -type d ! -name logs | sort)
PLEX="${PLEXES[$((SLURM_ARRAY_TASK_ID - 1))]}"
PNAME=$(basename "$PLEX")

echo "=== Task $SLURM_ARRAY_TASK_ID -> $PNAME ($(ls "$PLEX"/*.raw 2>/dev/null | wc -l) raw) ==="
for raw in "$PLEX"/*.raw; do
    [ -e "$raw" ] || continue
    base=$(basename "$raw" .raw)
    out="$PLEX/$base.mzML"
    [ -s "$out" ] && { echo "  SKIP $base"; continue; }
    echo "  CONVERT $base"
    "$TRFP" -i="$raw" -o="$PLEX" -f=2 -l=3      # f=2 plain mzML; no -g
    [ -s "$out" ] && echo "    OK $base" || echo "    FAIL $base"
done
echo "=== Done $PNAME ==="
```

Run:

```bash
sbatch msconvert.sh
```

### Verify (before deleting any raw)

```bash
find /scratch/maropakis.a/spectra -name '*.mzML.gz' | wc -l   # expect 0
find /scratch/maropakis.a/spectra -name '*.mzML'    | wc -l   # expect 302
find /scratch/maropakis.a/spectra -name '*.raw'     | wc -l   # must match mzML

# any raw with no matching mzML (re-convert list; should be empty):
comm -23 \
  <(find /scratch/maropakis.a/spectra -name '*.raw'  | sed 's/\.raw$//'  | sort) \
  <(find /scratch/maropakis.a/spectra -name '*.mzML' | sed 's/\.mzML$//' | sort)
```

Once counts match (302 = 302):

```bash
find /scratch/maropakis.a/spectra -name '*.raw' -delete
```

> If a prior run left gzipped output, normalize first — rename mislabeled
> plain files, gunzip the genuinely compressed ones:
> ```bash
> find /scratch/maropakis.a/spectra -name '*.mzML.gz' | while read -r f; do
>   file "$f" | grep -q 'gzip compressed' && gunzip "$f" || mv -- "$f" "${f%.gz}"
> done
> ```

---

## Step 2 — Run FragPipe (per plex)

`run_plexes.py` builds a per-plex `.workflow` + `.fp-manifest` and runs FragPipe
headless. It injects three keys per plex: `database.db-path`,
`fragger.fragger-path`, `tmtintegrator.channel_num`.

> **MSFragger is not in `tools/`** and is **not** found via `--config-tools-folder`.
> Its path must be written into each workflow via `fragger.fragger-path`. The
> `--msfragger-path` arg does this. (IonQuant only needed if
> `ionquant.run-ionquant=true` in a template — currently off; TMT-Integrator
> handles quant.)

`submit_fragpipe.sh`:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=fragpipe
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/maropakis.a/Frag_outputs/logs/fp_%A_%a.out
#SBATCH --error=/scratch/maropakis.a/Frag_outputs/logs/fp_%A_%a.err
set -euo pipefail
export JAVA_HOME=$HOME/bin/jdk-17.0.18+8
export PATH=$JAVA_HOME/bin:$PATH

PLEX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/maropakis.a/Frag_outputs/plex_list.txt)

python3 /home/maropakis.a/scripts/FragPipe/run_plexes.py \
  --spectra-root   /scratch/maropakis.a/spectra \
  --fasta-dir      /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
  --template-dir   /home/maropakis.a/scripts/FragPipe/templates \
  --out-dir        /scratch/maropakis.a/Frag_outputs \
  --fragpipe-bin   /home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe \
  --tools-folder   /home/maropakis.a/fragpipe/fragpipe-24.0/tools \
  --msfragger-path /home/maropakis.a/fragpipe/MSFragger-4.4.1/MSFragger-4.4.1.jar \
  --only "$PLEX" --run
```

Run a single plex first, then the full array:

```bash
sbatch --array=1-1  submit_fragpipe.sh    # test
sbatch --array=1-18 submit_fragpipe.sh    # full
```

> `plex_list.txt` must list exactly the 18 plex ids (one per line) and must not
> contain `logs`.

---

## Notes / gotchas

- `.mzML.gz` is not searchable by MSFragger — keep output plain (`-f=2`, no `-g`).
- Never delete `.raw` until the mzML/raw counts match.
- Stale `.cal`/`.index` files or a prior `combined/` dir can cause crashes — clear
  a plex's `results/<plex>/` before re-running it.
```
