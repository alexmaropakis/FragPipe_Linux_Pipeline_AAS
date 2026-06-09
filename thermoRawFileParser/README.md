# RAW → mzML File Conversion

Converts Thermo `.raw` to mzML for all TMT plexes on Explorer (MGHPCC), as a
SLURM array (one task per plex). Output is **plain** `.mzML` — gzipped mzML
isn't searchable by MSFragger, and skipping gzip avoids a separate
decompression step.

## Layout

```
$SPECTRA/<plex>/    # raw + mzML, one folder per plex
$SPECTRA/logs/      # conversion logs (NOT a plex)
```

## Prerequisite

[ThermoRawFileParser 2.0.0](https://github.com/compomics/ThermoRawFileParser) -- self-contained Linux build (no Mono required)

Reference: Hulstaert N, Shofstahl J, Sachsenberg T, Walzer M, Barsnes H, Martens L, Perez-Riverol Y: ThermoRawFileParser: Modular, Scalable, and Cross-Platform RAW File Conversion [PMID 31755270].

## Script — `msconvert.sh`

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
    "$TRFP" -i="$raw" -o="$PLEX" -f=2 -l=3      # f=2 plain mzML; no -g
    [ -s "$out" ] && echo "OK $base" || echo "FAIL $base"
done
```

Set `--array=1-N` to the number of plex folders (excluding `logs/`), then:

```bash
sbatch msconvert.sh
```

## Verify before deleting any raw

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

## Notes

- Keep mzML uncompressed (`-f=2`, no `-g`) — MSFragger can't read `.mzML.gz`.
- Never delete `.raw` until mzML/raw counts match.
- `logs/` sits under `$SPECTRA` but isn't a plex — exclude it from any
  directory-derived list, or array indexing is off by one and the last plex
  is skipped.
- If a prior run left gzipped/mislabeled output, normalize first:
  ```bash
  find $SPECTRA -name '*.mzML.gz' | while read -r f; do
    file "$f" | grep -q 'gzip compressed' && gunzip "$f" || mv -- "$f" "${f%.gz}"
  done
  ```
- mzML ≈ 1 GB/file — check `df -h /scratch` before bulk conversion.
```
