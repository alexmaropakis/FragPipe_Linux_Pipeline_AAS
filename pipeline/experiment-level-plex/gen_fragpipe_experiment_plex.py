#!/usr/bin/env python3
"""
One FragPipe run over several TMT plexes, bridged by TMT-Integrator (TMT10 to TMTpro18).

  RUN         = one workflow, manifest, FASTA and submit script.
  EXPERIMENT  = one plex: a name, its spectra dir, and a KEY selecting its sample_map rows
                (TMT plex, sample_ID, or box code e.g. TMT0276).

Per experiment: build annotation rows, convert spectra to .mzML, symlink into
<spectra-root>/<run>/<experiment>/.
Per run: set channel_num, pad annotations with <experiment>_EmptyN, write them to annotations/
and the staged dirs, then write manifest, workflow, and submit_<run>.sh.

Bridges (sample_name contains tmtintegrator.ref_tag) are renamed Bridge1, Bridge2, ... across
the run, so each is unique but shares the ref_tag prefix.

sample_map layout (one .xlsx per run):
  TMT plex | TMT channel | ParticipantID | Group | MQ | sample_name | sample_ID
     3     |    126      |    Bridge     | ...   | 1  |   Bridge     |   S3        <- experiment cortex_1
     4     |    126      |    Bridge     | ...   | 1  |   Bridge     |   S4        <- experiment cortex_2

Example:
  python gen_fragpipe_experiment_plex.py \\
    --run        cortex_tsumagari \\
    --species    mouse \\
    --workflow   /home/$USER/scripts/Search_gen/FragPipe/templates/TMT10_MS2_Val.workflow \\
    --sample-map /scratch/$USER/Dependencies/sample_map/sample_map_tsumagari_cortex.xlsx \\
    --experiment cortex_1 /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_1 3 \\
    --experiment cortex_2 /scratch/$USER/MQ_raw/Tsumagari_2023/cortex_2 4 \\
    --fasta      /scratch/$USER/Dependencies/FASTA_fragpipe/S9_cortex_tsumagari_fragpipe.fasta \\
    --out-dir    /scratch/$USER/Frag_outputs \\
    --spectra-root /scratch/$USER/spectra

Then: sbatch /scratch/$USER/Frag_outputs/submit/submit_cortex_tsumagari.sh

Note: make sure to "source activate  /projects/slavov/AM/envs/py39" before running this code

"""

import argparse
import glob
import os
import re
import shlex
import subprocess
import sys
from collections import Counter
import numpy as np
import pandas as pd
from psims.mzml import MzMLWriter
from pyteomics import mzxml

TMTPRO = ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C',
          '131N', '131C', '132N', '132C', '133N', '133C', '134N', '134C', '135N']
TMT10 = TMTPRO[:9] + ['131']
PLEX_CHANNELS = {10: TMT10, 11: TMTPRO[:11], 16: TMTPRO[:16], 18: TMTPRO}
ORD = {c: i for i, c in enumerate(TMT10 + TMTPRO[9:])}

SPECIES = ('human', 'mouse')
PLEX_SIZE_RE = re.compile(r'_(\d+)plex[._]', re.IGNORECASE)  # 'TMT0301_16plex_f001' -> 16
KEY_COLS = ('tmt_plex', 'sample_id', 'box')


def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def by_channel(rows):
    return sorted(rows, key=lambda r: ORD.get(r[0], 999))


def plex_size(files):
    """Plex size from the first '_<N>plex_' filename, else None."""
    for f in files:
        m = PLEX_SIZE_RE.search(os.path.basename(f))
        if m:
            return int(m.group(1))
    return None


# --- sample_map / annotation ---

def box_code(name):
    """'TMT0301_ID0131_16plex_02_127N' -> 'TMT0301'; 'TMT401_18plex' -> 'TMT0401'."""
    m = re.match(r'TMT0*(\d+)', str(name).strip()) if pd.notna(name) else None
    return f'TMT{int(m.group(1)):04d}' if m else None


def load_sample_map(path):
    df = pd.read_excel(path)
    df.columns = [re.sub(r'\s+', '_', str(c).strip().lower()) for c in df.columns]
    if not {'tmt_channel', 'sample_name'} <= set(df.columns):
        sys.exit(f'{path}: need tmt_channel + sample_name columns, got {list(df.columns)}')
    df = df.dropna(subset=['tmt_channel', 'sample_name'])
    for c in ('tmt_channel', 'sample_name'):
        df[c] = df[c].astype(str).str.strip()
    if 'unique_tmt_batch_name' in df.columns:  # several boxes can share one plex number
        df['box'] = df['unique_tmt_batch_name'].map(box_code)
    return df


def select_plex(df, key):
    """Rows whose tmt_plex, sample_id, or box matches KEY (case-insensitive)."""
    key = str(key).strip().lower()
    cols = [c for c in KEY_COLS if c in df.columns]
    for col in cols:
        sub = df[df[col].astype(str).str.strip().str.lower() == key]
        if len(sub):
            return sub
    sys.exit(f'sample_map: no rows for key {key!r} (searched: {cols})')


def annotation_rows(df, ref_tag, n):
    """(channel, name) rows; bridges -> <ref_tag><n>, other in-plex duplicates -> name_<channel>.
    Returns (rows, next n)."""
    rows = by_channel(zip(df['tmt_channel'], df['sample_name']))
    counts = Counter(name for _, name in rows)
    out = []
    for ch, name in rows:
        if ref_tag.lower() in name.lower():
            out.append((ch, f'{ref_tag}{n}'))
            n += 1
        else:
            out.append((ch, f'{name}_{ch}' if counts[name] > 1 else name))
    return out, n


def pad_rows(rows, size, experiment):
    """Fill unlabeled channels with <experiment>_EmptyN (TMT-Integrator needs size == channel_num)."""
    if len(rows) > size:
        print(f'  WARN: {experiment} has {len(rows)} channels > channel_num={size}; not padded')
        return rows
    used = {ch for ch, _ in rows}
    free = [c for c in PLEX_CHANNELS.get(size, TMTPRO[:size]) if c not in used][:size - len(rows)]
    return by_channel(rows + [(ch, f'{experiment}_Empty{i}') for i, ch in enumerate(free, 1)])


def write_annotation(rows, *paths):
    text = ''.join(f'{ch} {name}\n' for ch, name in rows)
    for p in paths:
        with open(p, 'w') as f:
            f.write(text)


# --- spectra ---

def mzxml_to_mzml(src, dst):
    """Pure-Python fallback when msconvert is missing. Writes <dst>.partial, renames on success.
    All header sections are required; without instrument_configuration_list MSFragger rejects it."""
    tmp = dst + '.partial'
    try:
        with mzxml.read(src, use_index=True) as reader, MzMLWriter(open(tmp, 'wb')) as w:
            w.controlled_vocabularies()
            w.file_description(file_contents=['MS1 spectrum', 'MS2 spectrum', 'centroid spectrum'],
                               source_files=[])
            sw = w.Software(id='pyteomics_psims_fallback', version='1.0',
                            params=['python-based mzXML to mzML converter'])
            w.software_list([sw])
            ic = w.InstrumentConfiguration(id='IC1', component_list=w.ComponentList([
                w.Source(order=1, params=['electrospray ionization']),
                w.Analyzer(order=2, params=['quadrupole']),
                w.Detector(order=3, params=['electron multiplier'])]))
            w.instrument_configuration_list([ic])
            w.data_processing_list([w.DataProcessing(
                [w.ProcessingMethod(order=1, software_reference=sw, params=['Conversion to mzML'])],
                id='DP1')])
            with w.run(id=os.path.splitext(os.path.basename(src))[0], instrument_configuration=ic), \
                    w.spectrum_list(count=len(reader)):
                for n, spec in enumerate(reader, 1):
                    mz = np.asarray(spec['m/z array'], dtype=np.float64)
                    inten = np.asarray(spec['intensity array'], dtype=np.float64)
                    level = int(spec.get('msLevel', 1))
                    prec = None
                    if level > 1 and spec.get('precursorMz'):
                        p = spec['precursorMz'][0]
                        prec = {'mz': float(p.get('precursorMz', 0.0)),
                                'intensity': float(p.get('precursorIntensity') or 0.0),
                                'activation': [p.get('activationMethod', 'HCD')]}
                        if p.get('precursorCharge'):
                            prec['charge'] = int(p['precursorCharge'])
                    w.write_spectrum(
                        mz, inten, id=f"scan={spec.get('num', n)}",
                        polarity='positive' if str(spec.get('polarity', '+')) == '+' else 'negative',
                        centroided=True, scan_start_time=float(spec.get('retentionTime', 0.0)),
                        precursor_information=prec,
                        params=[{'ms level': level}, {'total ion current': float(inten.sum())}])
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    os.rename(tmp, dst)


def collect_spectra(raw_dir, a):
    """One experiment's spectra as .mzML (plus .raw with --allow-raw); .mzXML is always converted."""
    ls = lambda ext: sorted(glob.glob(os.path.join(raw_dir, f'*.{ext}')))
    raws = ls('raw')
    if not (a.no_convert or a.allow_raw):
        for raw in raws:
            if not nonempty(os.path.splitext(raw)[0] + '.mzML'):
                # -f=2: plain indexed mzML (FragPipe skips .mzML.gz)
                subprocess.run([a.trfp, f'-i={raw}', f'-o={raw_dir}', '-f=2', '-l=3'], check=True)
                print(f'    converted {os.path.basename(raw)}')

    msconvert = shlex.split(a.msconvert)
    for f in ls('mzXML'):
        dst = os.path.splitext(f)[0] + '.mzML'
        if nonempty(dst):
            continue
        if msconvert:
            try:
                subprocess.run(msconvert + [f, '-o', raw_dir, '--mzML'], check=True)
            except FileNotFoundError:
                msconvert = None  # don't retry for the remaining files
        if not msconvert:
            mzxml_to_mzml(f, dst)
        print(f'    converted {os.path.basename(f)}' + ('' if msconvert else ' (Python fallback)'))
    return (raws if a.allow_raw else []) + ls('mzML')


def stage(files, dst):
    """Symlink FILES into DST, dropping stale .mzXML links; return the link paths."""
    os.makedirs(dst, exist_ok=True)
    for stale in glob.glob(os.path.join(dst, '*.mzXML')):
        os.remove(stale)
    links = []
    for f in files:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
        links.append(link)
    print(f'    staged {len(links)} spectra -> {dst}')
    return links


# --- run outputs ---

def patch_line(text, key, value):
    """Replace `key=...`, or append it if absent."""
    line = f'{key}={value}'
    out, n = re.subn(rf'^{re.escape(key)}=.*$', lambda _: line, text, flags=re.MULTILINE)
    return out if n else text.rstrip('\n') + f'\n{line}\n'


def read_ref_tag(text, default='Bridge'):
    m = re.search(r'^tmtintegrator\.ref_tag=(.+)$', text, re.MULTILINE)
    return m.group(1).strip() if m else default


SUBMIT_TEMPLATE = """\
#!/usr/bin/env bash
#SBATCH --job-name=fp_{run}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={ram}G
#SBATCH --time={time}
#SBATCH --output={logdir}/fp_{run}_%j.out
#SBATCH --error={logdir}/fp_{run}_%j.err
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
    ap = argparse.ArgumentParser(description='Prep one FragPipe run spanning multiple TMT plexes.')
    ap.add_argument('--run', required=True, help='run name, e.g. cortex_tsumagari')
    ap.add_argument('--species', required=True, choices=SPECIES)
    ap.add_argument('--experiment', required=True, action='append', nargs=3,
                    metavar=('NAME', 'RAW_DIR', 'PLEX'),
                    help='experiment name, its spectra dir, and its sample_map key '
                         '(TMT plex, sample_ID, or box code). Repeat per plex.')
    ap.add_argument('--sample-map', required=True, help='multi-plex sample_map .xlsx')
    ap.add_argument('--workflow', required=True, help='FragPipe .workflow template')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out-dir', required=True, help='Frag_outputs root')
    ap.add_argument('--spectra-root', required=True)
    ap.add_argument('--trfp', default=os.path.expanduser('~/thermoRawFileParser/ThermoRawFileParser'))
    ap.add_argument('--msconvert', default='msconvert',
                    help="msconvert command for .mzXML, e.g. 'singularity exec pwiz.sif msconvert'")
    ap.add_argument('--no-convert', action='store_true', help='skip .raw -> .mzML')
    ap.add_argument('--allow-raw', action='store_true', help='put .raw in the manifest unconverted')
    ap.add_argument('--fragpipe-bin', default=os.path.expandvars('/home/$USER/fragpipe/fragpipe-24.0/bin/fragpipe'))
    ap.add_argument('--tools-folder', default=os.path.expandvars('/home/$USER/fragpipe/fragpipe-24.0/tools'))
    ap.add_argument('--java-home', default=os.path.expanduser('~/bin/jdk-17.0.18+8'))
    ap.add_argument('--partition', default='short')
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--ram', type=int, default=64)
    ap.add_argument('--time', default='24:00:00')
    return ap.parse_args()


def main():
    a = parse_args()
    run = a.run.strip().lower()
    if not os.path.isfile(a.fasta):
        sys.exit(f'--fasta not found: {a.fasta}')
    out = {d: os.path.join(a.out_dir, d) for d in ('workflows', 'manifests', 'annotations', 'submit', 'logs')}
    workdir = os.path.abspath(os.path.join(a.out_dir, 'results', run))
    for d in [*out.values(), workdir]:
        os.makedirs(d, exist_ok=True)

    template = open(a.workflow).read()
    ref_tag = read_ref_tag(template)
    smap = load_sample_map(a.sample_map)
    print(f'[{run}] species={a.species}  {len(a.experiment)} experiments  bridges -> {ref_tag}N')

    exps, entries, next_bridge = {}, [], 1  # entries: (staged path, experiment) for the manifest
    for name, raw_dir, key in a.experiment:
        name = name.strip()
        if name in exps:
            sys.exit(f'duplicate experiment name {name!r}')
        if not os.path.isdir(raw_dir):
            sys.exit(f'[{name}] raw_dir not found: {raw_dir}')
        print(f'  experiment {name!r}  key={key}  raw_dir={raw_dir}')

        first = next_bridge
        rows, next_bridge = annotation_rows(select_plex(smap, key), ref_tag, first)
        bridges = ', '.join(f'{ref_tag}{i}' for i in range(first, next_bridge))
        print(f'    {len(rows)} channels, bridges: {bridges or "NONE (cannot be bridged)"}')

        files = collect_spectra(raw_dir, a)
        if not files:
            sys.exit(f'[{name}] no .mzML/.mzXML/.raw found')
        size = plex_size(files)
        if size is None:
            print(f'    WARN: no _<N>plex_ in filenames; using {len(rows)} labeled channels')

        staged_dir = os.path.join(a.spectra_root, run, name)
        entries += [(p, name) for p in stage(files, staged_dir)]
        exps[name] = dict(rows=rows, size=size, dir=staged_dir)

    # one workflow -> one channel_num; filename plex size beats labeled count
    sizes = {n: e['size'] for n, e in exps.items() if e['size']}
    if len(set(sizes.values())) > 1:
        print(f'  WARN: experiments differ in plex size {sizes}; split the run by plex size')
    channels = max(e['size'] or len(e['rows']) for e in exps.values())

    for name, e in exps.items():
        write_annotation(pad_rows(e['rows'], channels, name),
                         os.path.join(out['annotations'], f'{run}_{name}_annotation.txt'),
                         os.path.join(e['dir'], 'annotation.txt'))  # TMT-Integrator reads this one

    wf_path = os.path.join(out['workflows'], f'{run}.workflow')
    mf_path = os.path.join(out['manifests'], f'{run}.fp-manifest')
    tbl_path = os.path.join(out['annotations'], f'{run}_samples.txt')
    submit_path = os.path.join(out['submit'], f'submit_{run}.sh')

    wf = patch_line(template, 'database.db-path', os.path.abspath(a.fasta))
    wf = patch_line(wf, 'tmtintegrator.channel_num', str(channels))
    wf = patch_line(wf, 'tmtintegrator.output', os.path.join(workdir, 'tmt-report'))  # must not be blank
    with open(wf_path, 'w') as f:
        f.write(wf)
    with open(mf_path, 'w') as f:
        f.writelines(f'{os.path.abspath(p)}\t{n}\t1\tDDA\n' for p, n in entries)
    with open(tbl_path, 'w') as f:
        f.write('file\texperiment\n')
        f.writelines(f'{os.path.basename(p)}\t{n}\n' for p, n in entries)
    with open(submit_path, 'w') as f:
        f.write(SUBMIT_TEMPLATE.format(
            run=run, partition=a.partition, threads=a.threads, ram=a.ram, time=a.time,
            logdir=out['logs'], java_home=a.java_home, fragpipe_bin=a.fragpipe_bin,
            workflow=wf_path, manifest=mf_path, workdir=workdir, tools_folder=a.tools_folder))
    os.chmod(submit_path, 0o755)

    print(f'  channel_num = {channels}\n  workflow -> {wf_path}\n  manifest -> {mf_path} '
          f'({len(entries)} files, {len(exps)} experiments)\n  samples  -> {tbl_path}'
          f'\n  submit   -> {submit_path}\n\nNext: sbatch {submit_path}')


if __name__ == '__main__':
    main()
