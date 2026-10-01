#!/usr/bin/env python3
"""
Prepare one TMT plex for FragPipe (TMT-labelled data only).

  1. annotation.txt from the plex's sample_map, sorted by channel (in-plex duplicates -> name_<channel>)
  2. .raw -> .mzML with ThermoRawFileParser, skipping existing (--no-convert to skip)
  3. <spectra-root>/<plex>/ with symlinked .mzML + annotation.txt
  4. workflow (db-path, channel_num) and .fp-manifest
  5. submit_<plex>.sh

Example:
  python gen_fragpipe_plex.py /scratch/maropakis.a/MQ_raw/Ping_2018/ACG/b1 \\
    --plex       acgb1 \\
    --species    human \\
    --channels   10 \\
    --workflow   /home/maropakis.a/scripts/FragPipe/templates/TMT10_MS3_Val.workflow \\
    --sample-map /scratch/maropakis.a/Dependencies/sample_map/acgb1.xlsx \\
    --fasta-dir  /scratch/maropakis.a/Dependencies/FASTA_fragpipe \\
    --out-dir    /scratch/maropakis.a/Frag_outputs \\
    --spectra-root /scratch/maropakis.a/spectra

Then: sbatch /scratch/maropakis.a/Frag_outputs/submit/submit_acgb1.sh
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import Counter

import pandas as pd

CHANNEL_ORDER = ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C',
                 '131', '131N', '131C', '132N', '132C', '133N', '133C', '134N', '134C', '135N']
ORD = {c: i for i, c in enumerate(CHANNEL_ORDER)}
SPECIES = ('human', 'mouse')


def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def write_annotation(sample_map, out_path, expected=None):
    """sample_map .xlsx -> annotation.txt; returns the channel count."""
    df = pd.read_excel(sample_map)
    df.columns = [re.sub(r'\s+', '_', str(c).strip().lower()) for c in df.columns]
    if not {'tmt_channel', 'sample_name'} <= set(df.columns):
        sys.exit(f'{sample_map}: need tmt_channel + sample_name columns, got {list(df.columns)}')
    df = df.dropna(subset=['tmt_channel', 'sample_name'])
    chans = df['tmt_channel'].astype(str).str.strip()
    names = df['sample_name'].astype(str).str.strip()
    counts = Counter(names)
    rows = sorted(((ch, f'{n}_{ch}' if counts[n] > 1 else n) for ch, n in zip(chans, names)),
                  key=lambda r: ORD.get(r[0], 999))
    with open(out_path, 'w') as f:
        f.writelines(f'{ch} {n}\n' for ch, n in rows)
    if expected is not None and len(rows) != expected:
        print(f'  WARN: annotation has {len(rows)} channels but --channels {expected}')
    return len(rows)


def convert_raws(raw_dir, trfp):
    """Convert each .raw to .mzML in place (skip existing); returns .mzML paths."""
    raws = sorted(glob.glob(os.path.join(raw_dir, '*.raw')))
    if not raws:
        mzmls = sorted(glob.glob(os.path.join(raw_dir, '*.mzML')))
        print(f'  no .raw in {raw_dir}; found {len(mzmls)} existing .mzML')
        return mzmls
    mzmls = []
    for raw in raws:
        out = os.path.splitext(raw)[0] + '.mzML'
        base = os.path.basename(os.path.splitext(raw)[0])
        if nonempty(out):
            print(f'  SKIP convert {base} (mzML exists)')
        else:
            # -f=2: plain indexed mzML (FragPipe skips .mzML.gz)
            subprocess.run([trfp, f'-i={raw}', f'-o={raw_dir}', '-f=2', '-l=3'], check=True)
            print(f'  converted {base}')
        if os.path.exists(out):
            mzmls.append(out)
    return mzmls


def find_fasta(fasta_dir, plex):
    """S1_ACGB1_fragpipe.fasta / S9_cortex_keele_MTP_fragpipe.fasta -> match plex (case-insensitive)."""
    for p in sorted(glob.glob(os.path.join(fasta_dir, '*_fragpipe.fasta'))):
        label = re.sub(r'^S\d+_', '', os.path.basename(p), flags=re.I)
        if re.sub(r'(?:_MTP)?_fragpipe\.fasta$', '', label, flags=re.I).lower() == plex:
            return p
    sys.exit(f'no *_fragpipe.fasta in {fasta_dir} for plex {plex!r}')


def patch_line(text, key, value):
    """Replace `key=...`, or append it if absent."""
    line = f'{key}={value}'
    out, n = re.subn(rf'^{re.escape(key)}=.*$', lambda _: line, text, flags=re.MULTILINE)
    return out if n else text.rstrip('\n') + f'\n{line}\n'


SUBMIT_TEMPLATE = """\
#!/usr/bin/env bash
#SBATCH --job-name=fp_{plex}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={ram}G
#SBATCH --time={time}
#SBATCH --output={logdir}/fp_{plex}_%j.out
#SBATCH --error={logdir}/fp_{plex}_%j.err
set -euo pipefail
export JAVA_HOME={java_home}
export PATH=$JAVA_HOME/bin:$PATH

{fragpipe_bin} --headless \\
  --workflow {workflow} \\
  --manifest {manifest} \\
  --workdir  {workdir} \\
  --threads  {threads} \\
  --ram      {ram} \\
  --config-tools-folder {tools_folder}
"""


def parse_args():
    ap = argparse.ArgumentParser(description='Prep one FragPipe plex end-to-end + submit script.')
    ap.add_argument('raw_dir', help="dir holding this plex's .raw (or pre-made .mzML)")
    ap.add_argument('--plex', required=True, help='plex token, e.g. acgb1 / pooled / aorta')
    ap.add_argument('--species', required=True, choices=SPECIES)
    ap.add_argument('--channels', type=int, default=None,
                    help='expected channel count; only checked (count comes from the sample_map)')
    ap.add_argument('--workflow', required=True, help='FragPipe .workflow template')
    ap.add_argument('--sample-map', required=True, help="this plex's sample_map .xlsx")
    ap.add_argument('--fasta-dir', required=True, help='dir of *_fragpipe.fasta')
    ap.add_argument('--out-dir', required=True, help='Frag_outputs root')
    ap.add_argument('--spectra-root', required=True)
    ap.add_argument('--trfp', default=os.path.expanduser('~/thermoRawFileParser/ThermoRawFileParser'))
    ap.add_argument('--no-convert', action='store_true', help='skip .raw -> .mzML')
    ap.add_argument('--fragpipe-bin', default='/home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe')
    ap.add_argument('--tools-folder', default='/home/maropakis.a/fragpipe/fragpipe-24.0/tools')
    ap.add_argument('--java-home', default=os.path.expanduser('~/bin/jdk-17.0.18+8'))
    ap.add_argument('--partition', default='short')
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--ram', type=int, default=64)
    ap.add_argument('--time', default='24:00:00')
    return ap.parse_args()


def main():
    a = parse_args()
    plex = a.plex.strip().lower()
    if not os.path.isdir(a.raw_dir):
        sys.exit(f'raw_dir not found: {a.raw_dir}')
    out = {d: os.path.join(a.out_dir, d) for d in ('workflows', 'manifests', 'annotations', 'submit')}
    for d in out.values():
        os.makedirs(d, exist_ok=True)
    print(f'[{plex}] species={a.species} channels={a.channels or "auto"}')

    annot_path = os.path.join(out['annotations'], f'{plex}_annotation.txt')
    channels = write_annotation(a.sample_map, annot_path, a.channels)
    print(f'  {channels} channels')

    if a.no_convert:
        mzmls = sorted(glob.glob(os.path.join(a.raw_dir, '*.mzML')))
        print(f'  --no-convert: using {len(mzmls)} existing .mzML')
    else:
        mzmls = convert_raws(a.raw_dir, a.trfp)
    if not mzmls:
        sys.exit(f'{plex}: no mzML to stage')

    stage_dir = os.path.join(a.spectra_root, plex)
    os.makedirs(stage_dir, exist_ok=True)
    for f in mzmls:
        link = os.path.join(stage_dir, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
    shutil.copy(annot_path, os.path.join(stage_dir, 'annotation.txt'))
    print(f'  staged {len(mzmls)} mzML + annotation -> {stage_dir}')
    staged = sorted(glob.glob(os.path.join(stage_dir, '*.mzML')))

    fasta = find_fasta(a.fasta_dir, plex)
    wf_path = os.path.join(out['workflows'], f'{plex}.workflow')
    mf_path = os.path.join(out['manifests'], f'{plex}.fp-manifest')
    submit_path = os.path.join(out['submit'], f'submit_{plex}.sh')
    workdir = os.path.join(a.out_dir, 'results', plex)
    logdir = os.path.join(a.out_dir, 'logs')
    for d in (workdir, logdir):
        os.makedirs(d, exist_ok=True)

    wf = patch_line(open(a.workflow).read(), 'database.db-path', os.path.abspath(fasta))
    wf = patch_line(wf, 'tmtintegrator.channel_num', str(channels))
    with open(wf_path, 'w') as f:
        f.write(wf)
    with open(mf_path, 'w') as f:
        f.writelines(f'{os.path.abspath(s)}\t{plex}\t1\tDDA\n' for s in staged)
    with open(submit_path, 'w') as f:
        f.write(SUBMIT_TEMPLATE.format(
            plex=plex, partition=a.partition, threads=a.threads, ram=a.ram, time=a.time,
            logdir=logdir, java_home=a.java_home, fragpipe_bin=a.fragpipe_bin,
            workflow=wf_path, manifest=mf_path, workdir=workdir, tools_folder=a.tools_folder))
    os.chmod(submit_path, 0o755)

    print(f'  workflow -> {wf_path}\n  manifest -> {mf_path} ({len(staged)} files)'
          f'\n  fasta = {os.path.basename(fasta)}\n  submit   -> {submit_path}\n\nNext: sbatch {submit_path}')


if __name__ == '__main__':
    main()
