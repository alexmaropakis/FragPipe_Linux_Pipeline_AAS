#!/usr/bin/env python3

"""
Script to run FragPipe at the experiment level with TMT11plex Bridging. 
Use this when several TMT plexes must be searched/quantified together. 

  RUN         = the whole FragPipe job (one workflow, one manifest, one submit script, one FASTA).
  EXPERIMENT  = one TMT plex inside that run. You name it (cortex_1, cortex_2, hippocampus-1 ...),
                point it at the folder holding that plex's .raw/.mzML/.mzXML, and give the plex KEY
                that selects that plex's rows from the shared sample_map (its 'TMT plex' or 'sample_ID').

For each experiment:
1. Builds annotation.txt from sample_map for that plex
2. Converts/stages spectra into <spectra-root>/<run>/<experiment>/ (symlinked mzML/raw + that experiment's annotation.txt)
3. Writes manifest listing every file across all experiments
   -- ONLY .mzML files (or .raw, if --allow-raw is passed) ever go into the manifest / staged
   dirs. .mzXML is always converted away and never written out, even if a stale .mzXML symlink
   is sitting in an old staged directory from a previous run.
4. Patches workflow template 
5. Writes submit_<run>.sh 

Multi-plex sample_map layout (one .xlsx for the whole run), e.g. sample_map_tsumagari_cortex.xlsx:
  TMT plex | TMT channel | ParticipantID | Group | MQ | sample_name | sample_ID
     3     |    126      |    Bridge     | ...   | 1  |   Bridge     |   S3        <- experiment cortex_1
     4     |    126      |    Bridge     | ...   | 1  |   Bridge     |   S4        <- experiment cortex_2

Example usage: (cortex_1 = plex 3, cortex_2 = plex 4)
  python gen_fragpipe_experiment_plex.py \
    --run        cortex_tsumagari \
    --species    mouse \
    --workflow   /home/maropakis.a/scripts/Search_gen/FragPipe/templates/TMT10_MS2_Val.workflow \
    --sample-map /scratch/maropakis.a/Dependencies/sample_map/sample_map_tsumagari_cortex.xlsx \
    --experiment cortex_1 /scratch/maropakis.a/MQ_raw/Tsumagari_2023/cortex_1 3 \
    --experiment cortex_2 /scratch/maropakis.a/MQ_raw/Tsumagari_2023/cortex_2 4 \
    --fasta      /scratch/maropakis.a/Dependencies/FASTA_fragpipe/S9_cortex_tsumagari_fragpipe.fasta \
    --out-dir    /scratch/maropakis.a/Frag_outputs \
    --spectra-root /scratch/maropakis.a/spectra

Then: sbatch /scratch/maropakis.a/Frag_outputs/submit/submit_cortex_tsumagari.sh               

"""

import argparse
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter

# Header parsing 
# canonical TMTpro/TMT channel order (covers TMT10 through TMTpro18)
CHANNEL_ORDER = ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C',
                 '131', '131N', '131C', '132N', '132C', '133N', '133C', '134N', '134C', '135N']
ORD = {c: i for i, c in enumerate(CHANNEL_ORDER)}

SPECIES_TAG = {
    'human': ('Homo sapiens', 9606),
    'mouse': ('Mus musculus', 10090),
}

# spectra extensions this script will look for on disk, besides converting .raw. NOTE: this is
# only used to *find* candidate spectra to feed into the convert/ensure_mzml pipeline -- it is
# NOT the set of extensions allowed into the final manifest (see ALLOWED_MANIFEST_EXTS below).
SPECTRA_EXTS = ('mzML', 'mzXML')

# 'TMT0301_16plex_f001.mzXML' -> 16. This is the true physical TMT reagent size (channel_num),
# which is NOT the same as how many of those channels happen to be labeled in a given cohort's
# slice of a shared/bridged sample_map (see detect_plex_size below).
PLEX_SIZE_RE = re.compile(r'_(\d+)plex[._]', re.IGNORECASE)

def detect_plex_size(files):
    """Parse the TMT reagent channel count from a spectra filename. Returns None if none of
    the given files matches the '_<N>plex_' naming convention."""
    for f in files:
        m = PLEX_SIZE_RE.search(os.path.basename(f))
        if m:
            return int(m.group(1))
    return None

# Helper functions
def norm_col(col):
    """'Sample name'/'TMT channel'/'sample_ID' -> 'sample_name'/'tmt_channel'/'sample_id'."""
    return re.sub(r'\s+', '_', str(col).strip().lower())

def glob_spectra(dir_path):
    """All .mzML/.mzXML files in dir_path, sorted."""
    files = []
    for ext in SPECTRA_EXTS:
        files.extend(glob.glob(os.path.join(dir_path, f'*.{ext}')))
    return sorted(files)

def box_from_batch_name(name):
    """'TMT0301_ID0131_16plex_02_127N' -> 'TMT0301'; 'TMT401_18plex_127N' -> 'TMT0401'
    (zero-padded to 4 digits, matching the actual mzXML filenames)."""
    import pandas as pd
    if pd.isna(name):
        return None
    m = re.match(r'TMT0*(\d+)', str(name).strip())
    return f'TMT{int(m.group(1)):04d}' if m else None

def load_sample_map(path):
    """Read the (possibly multi-plex) sample_map .xlsx once; normalize columns + string cells."""
    import pandas as pd
    df = pd.read_excel(path)
    df.columns = [norm_col(c) for c in df.columns]
    if not {'tmt_channel', 'sample_name'} <= set(df.columns):
        sys.exit(f'{path}: need tmt_channel + sample_name columns, got {list(df.columns)}')
    df = df.dropna(subset=['tmt_channel', 'sample_name'])
    df['tmt_channel'] = df['tmt_channel'].astype(str).str.strip()
    df['sample_name'] = df['sample_name'].astype(str).str.strip()
    if 'unique_tmt_batch_name' in df.columns:
        # derived physical-box code, e.g. 'TMT0276' -- lets select_plex isolate one box's own
        # channel->sample rows even when several boxes share one 'TMT plex' number (matched_tissues)
        df['box'] = df['unique_tmt_batch_name'].map(box_from_batch_name)
    return df

def select_plex(df, key):
    """Return the sample_map rows for one plex/box, matching KEY against tmt_plex, sample_id, or box.

    KEY is matched case-insensitively against whichever column has it, so '3', 'S3', and 'TMT0276'
    are all valid. Prefer passing a box code (e.g. 'TMT0276') over a plex number whenever a sample_map
    can have several distinct physical boxes sharing one 'TMT plex' value (matched_tissues) -- a plex
    number there is not enough to pick out one consistent channel->sample mapping. A run with a
    single-plex sample_map (no plex column) may pass key=None to take the whole table.
    """
    if key is None:
        return df
    key = str(key).strip().lower()
    for col in ('tmt_plex', 'sample_id', 'box'):
        if col in df.columns:
            sub = df[df[col].astype(str).str.strip().str.lower() == key]
            if len(sub):
                return sub
    have = [c for c in ('tmt_plex', 'sample_id', 'box') if c in df.columns] or ['<none>']
    sys.exit(f'sample_map: no rows for plex/box key {key!r} (searched columns: {", ".join(have)})')

def write_annotation(df, out_path, expected=None):
    """One plex's rows -> FragPipe annotation.txt ('<channel> <sample>' per line).
    Returns (row_count, {physical channels used}) -- the latter lets pad_annotation() later
    fill in any remaining channels of the plex with Empty placeholders.

    Repeated labels within the plex (e.g. a Bridge on two channels) are disambiguated with _<ch>,
    matching Alex_gen_fragpipe.py. Rows are ordered by canonical TMT channel.
    """
    counts = Counter(df['sample_name'])
    rows = [(ch, f'{name}_{ch}' if counts[name] > 1 else name)
            for ch, name in zip(df['tmt_channel'], df['sample_name'])]
    rows.sort(key=lambda x: ORD.get(x[0], 999))
    with open(out_path, 'w') as f:
        for ch, name in rows:
            f.write(f'{ch} {name}\n')
    if expected is not None and len(rows) != expected:
        print(f'  WARN: {os.path.basename(out_path)} has {len(rows)} channels but expected {expected}')
    return len(rows), {ch for ch, _ in rows}

def pad_annotation(annot_path, used_channels, target_count, experiment):
    """Ensure annot_path ends up with exactly TARGET_COUNT channel lines.

    TMTIntegrator requires the annotation file's channel count to exactly equal channel_num on
    the Quant (Isobaric) tab -- a mismatch throws 'Number of the samples in the annotation file
    does not match the number of channels...'. This happens whenever the sample_map only labels
    a subset of a plex's physical channels (e.g. 16 of 18 used), since write_annotation() only
    ever writes rows for samples that are actually present.

    Fills any unlabeled physical channels (taken from the first TARGET_COUNT slots of the
    canonical CHANNEL_ORDER -- i.e. that plex's real channel set) with placeholder rows named
    '<experiment>_Empty1', '<experiment>_Empty2', ... . The experiment prefix is required, not
    cosmetic: TMTIntegrator also demands sample names be unique across EVERY experiment's
    annotation file in the whole run, so plain 'Empty1'/'Empty2' would collide the moment more
    than one experiment needs padding. EmptyN numbering is always contiguous starting at 1 within
    one annotation file, in canonical channel order.

    The full file (existing rows + new Empty rows) is re-sorted by canonical channel order and
    rewritten, rather than appending Empty rows at the end -- otherwise they'd land after
    whatever channel happened to be labeled last (e.g. '126 ... Empty1' printed after '134N ...'),
    which reads as if the file were out of order even though every channel is still present once.
    """
    have = len(used_channels)
    if have >= target_count:
        if have > target_count:
            print(f'  WARN: {os.path.basename(annot_path)} has {have} channels, more than '
                  f'channel_num={target_count}; leaving as-is -- check this manually')
        return
    plex_channels = CHANNEL_ORDER[:target_count]
    unused = [c for c in plex_channels if c not in used_channels]
    n_missing = target_count - have
    fill = unused[:n_missing]
    if len(fill) < n_missing:
        # shouldn't normally happen (would mean target_count > len(CHANNEL_ORDER)), but don't
        # silently under-fill -- fall back to synthetic channel labels so the count still matches
        fill += [f'Channel{i}' for i in range(len(fill) + 1, n_missing + 1)]

    with open(annot_path) as f:
        existing = [tuple(line.split(maxsplit=1)) for line in f if line.strip()]
    existing = [(ch, name.rstrip('\n')) for ch, name in existing]

    new_rows = [(ch, f'{experiment}_Empty{i}') for i, ch in enumerate(fill, start=1)]
    all_rows = sorted(existing + new_rows, key=lambda r: ORD.get(r[0], 999))

    with open(annot_path, 'w') as f:
        for ch, name in all_rows:
            f.write(f'{ch} {name}\n')
    print(f'  padded {os.path.basename(annot_path)}: +{len(fill)} Empty channel(s) '
          f'({have} -> {target_count}) to match channel_num, re-sorted by channel order')


def convert_raws(raw_dir, trfp, allow_raw=False):
    """Convert every .raw in raw_dir to plain mzML in place; skip existing. Return mzML/mzXML/raw paths.

    If allow_raw is True, .raw files found are NOT converted -- they're returned as-is and will
    be staged/manifested directly (TMTIntegrator can read .raw natively). Any .mzML/.mzXML already
    present alongside them is still picked up too (mixed dirs are handled downstream by ensure_mzml).
    """
    raws = sorted(glob.glob(os.path.join(raw_dir, '*.raw')))
    if not raws:
        # why: raw_dir may already hold mzML/mzXML (pre-converted, or e.g. sorted by plex from an
        # archive that only ships mzXML); fall through to collect those instead of converting.
        mzmls = glob_spectra(raw_dir)
        print(f'    no .raw in {raw_dir}; found {len(mzmls)} existing .mzML/.mzXML')
        return mzmls
    if allow_raw:
        print(f'    --allow-raw: keeping {len(raws)} .raw file(s) uncoverted in {raw_dir}')
        return raws + glob_spectra(raw_dir)
    mzmls = []
    for raw in raws:
        base = os.path.splitext(os.path.basename(raw))[0]
        out = os.path.join(raw_dir, f'{base}.mzML')
        if os.path.getsize(out) if os.path.exists(out) else 0:
            print(f'    SKIP convert {base} (mzML exists)')
        else:
            # why: -f=2 = plain indexed mzML; -g would gzip and FragPipe silently skips .mzML.gz.
            subprocess.run([trfp, f'-i={raw}', f'-o={raw_dir}', '-f=2', '-l=3'], check=True)
            print(f'    converted {base}')
        if os.path.exists(out):
            mzmls.append(out)
    return mzmls

def convert_mzxml_python(mzxml_path, dst):
    """Pure-Python .mzXML -> .mzML fallback (no ProteoWizard/msconvert needed), via pyteomics+psims.

    Writes to a temp '<dst>.partial' path and only os.rename()'s it to DST once the conversion is
    fully complete -- so a killed/interrupted run never leaves a truncated file sitting at DST,
    where callers' exists()+size skip-checks would otherwise mistake it for a finished conversion.
    Produces a functionally complete mzML (spectra, MS level, RT, polarity, precursor mz/charge/
    intensity) with the minimal but complete header section chain mzML actually requires.

    psims's MzMLWriter enforces mzML's section ordering as a state machine: controlled_vocabularies
    -> file_description -> software_list -> instrument_configuration_list -> data_processing_list
    -> run. Skipping a section (as an earlier version of this function did for
    instrument_configuration_list) doesn't just warn -- it leaves <run> without a valid
    defaultInstrumentConfigurationRef, which is REQUIRED by the mzML schema and is exactly what
    made FragPipe's CheckCentroid step (and MSFragger itself) refuse to open the file with
    'RunHeaderParsingException: Could not find "defaultInstrumentConfigurationRef" attribute'.
    file_description/software_list/instrument_configuration_list/data_processing_list are
    therefore all written every time, in this exact order, even though their *content* is
    intentionally generic/minimal (mzXML's own metadata doesn't reliably expose real instrument
    or software details, so this doesn't assert specifics we don't actually know). None of this
    affects identification or quantification, which don't depend on this header metadata.

    Lazily imports numpy/pyteomics/psims -- ONLY when this fallback actually runs. This matters:
    pyteomics pulls in psims, which pulls in h5py (via its mzmlb writer submodule), which has a
    history of NumPy 1.x/2.x ABI mismatches against whatever numpy happens to be on PATH in a
    given conda env. Since ensure_mzml() only calls this function when msconvert itself can't be
    found, a working --msconvert means this whole fragile import chain is never touched at all.
    """
    import numpy as np
    from pyteomics import mzxml
    from psims.mzml import MzMLWriter

    tmp = dst + '.partial'
    with mzxml.read(mzxml_path, use_index=True) as reader:
        try:
            n_total = len(reader)
        except TypeError:
            n_total = None
        try:
            with MzMLWriter(open(tmp, 'wb')) as writer:
                writer.controlled_vocabularies()

                writer.file_description(
                    file_contents=['MS1 spectrum', 'MS2 spectrum', 'centroid spectrum'],
                    source_files=[],
                )

                sw = writer.Software(id='pyteomics_psims_fallback', version='1.0',
                                      params=['python-based mzXML to mzML converter'])
                writer.software_list([sw])

                ic = writer.InstrumentConfiguration(
                    id='IC1',
                    component_list=writer.ComponentList([
                        writer.Source(order=1, params=['electrospray ionization']),
                        writer.Analyzer(order=2, params=['quadrupole']),
                        writer.Detector(order=3, params=['electron multiplier']),
                    ]),
                )
                writer.instrument_configuration_list([ic])

                dp = writer.DataProcessing(
                    [writer.ProcessingMethod(order=1, software_reference=sw,
                                              params=['Conversion to mzML'])],
                    id='DP1',
                )
                writer.data_processing_list([dp])

                with writer.run(id=os.path.splitext(os.path.basename(mzxml_path))[0],
                                 instrument_configuration=ic):
                    with writer.spectrum_list(count=n_total):
                        for n, spec in enumerate(reader, start=1):
                            mz = np.asarray(spec['m/z array'], dtype=np.float64)
                            inten = np.asarray(spec['intensity array'], dtype=np.float64)
                            ms_level = int(spec.get('msLevel', 1))
                            rt = float(spec.get('retentionTime', 0.0))  # pyteomics gives minutes
                            scan_num = spec.get('num', n)
                            polarity = 'positive' if str(spec.get('polarity', '+')) == '+' else 'negative'

                            precursor_info = None
                            plist = spec.get('precursorMz')
                            if ms_level > 1 and plist:
                                p = plist[0]
                                precursor_info = {
                                    'mz': float(p.get('precursorMz', 0.0)),
                                    'intensity': float(p.get('precursorIntensity', 0.0) or 0.0),
                                    'activation': [p.get('activationMethod', 'HCD')],
                                }
                                if p.get('precursorCharge'):
                                    precursor_info['charge'] = int(p['precursorCharge'])

                            writer.write_spectrum(
                                mz, inten,
                                id=f'scan={scan_num}',
                                polarity=polarity,
                                centroided=True,
                                scan_start_time=rt,
                                precursor_information=precursor_info,
                                params=[{'ms level': ms_level},
                                        {'total ion current': float(inten.sum()) if inten.size else 0.0}],
                            )
        except BaseException:
            # why: SIGKILL itself can't be caught, but this covers Ctrl-C / real errors; the
            # tmp-name + rename-on-success pattern below is what actually protects against SIGKILL.
            if os.path.exists(tmp):
                os.remove(tmp)
            raise
    os.rename(tmp, dst)  # only appears at the real name once fully written

def ensure_mzml(files, msconvert_cmd):
    """TMTIntegrator only accepts .mzML/.raw/.d -- NOT .mzXML. Convert any .mzXML in FILES to
    .mzML in place (skip if already converted) and return the list with .mzXML paths swapped
    for their .mzML counterpart. .raw/.d/.mzML pass through untouched.

    Tries `msconvert_cmd` (ProteoWizard msconvert, possibly wrapped e.g. in singularity/apptainer/
    wine) first -- it's the standard, most complete converter. If that command can't be found at
    all (no ProteoWizard install, no container runtime available), automatically falls back to a
    pure-Python conversion via pyteomics+psims for the rest of the run instead of failing outright.
    """
    cmd_prefix = shlex.split(msconvert_cmd)
    use_python_fallback = False
    out = []
    for f in files:
        if not f.lower().endswith('.mzxml'):
            out.append(f)
            continue
        d = os.path.dirname(f)
        base = os.path.splitext(os.path.basename(f))[0]
        dst = os.path.join(d, f'{base}.mzML')
        if os.path.exists(dst) and os.path.getsize(dst):
            print(f'    SKIP convert {base} (.mzML exists)')
            out.append(dst)
            continue

        if not use_python_fallback:
            try:
                subprocess.run(cmd_prefix + [f, '-o', d, '--mzML'], check=True)
                print(f'    converted {base}.mzXML -> .mzML (msconvert)')
                out.append(dst)
                continue
            except FileNotFoundError:
                use_python_fallback = True
                print(f"    msconvert command not found ({msconvert_cmd!r}) -- falling back to "
                      f'built-in Python mzXML->mzML conversion (pyteomics+psims) for this and any '
                      f'remaining .mzXML files')

        convert_mzxml_python(f, dst)
        print(f'    converted {base}.mzXML -> .mzML (python fallback)')
        out.append(dst)
    return out

def stage_experiment(mzmls, spectra_root, run, experiment):
    """Symlink one experiment's files (already resolved to .mzML/.raw by this point -- never
    .mzXML) into <spectra-root>/<run>/<experiment>/.

    NOTE: this does NOT copy annotation.txt in -- that must happen via copy_annotation_to_stage()
    *after* pad_annotation() has run (see main()). Copying it here, before channel_num is even
    known for the whole run, was a real bug: TmtIntegrator reads annotation.txt from THIS staged
    directory, not from the annotations/ dir, so an early copy meant the *unpadded* file was the
    one actually used, no matter how correctly pad_annotation() patched the other copy later.

    Returns (staged_paths, dst) -- the exact staged symlink paths, and the staged dir itself, so
    callers don't need to re-glob the directory. Re-globbing was a separate real bug: a staged
    dir left over from an older run/version of this script could still hold stale .mzXML symlinks
    (from before ensure_mzml swapped them to .mzML), and a directory glob for *.mzXML would
    silently pick those back up into the manifest even though the CURRENT run never asked for them.

    As a belt-and-suspenders cleanup, any stale *.mzXML symlink already sitting in DST that
    isn't one of MZMLS's own source files is removed here too.
    """
    dst = os.path.join(spectra_root, run, experiment)
    os.makedirs(dst, exist_ok=True)

    staged_paths = []
    for f in mzmls:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
        staged_paths.append(link)

    # belt-and-suspenders: purge any leftover .mzXML symlinks in this staged dir from a prior
    # run/version so they can never leak into this or any future manifest via a directory glob.
    for stale in glob.glob(os.path.join(dst, '*.mzXML')):
        if os.path.realpath(stale) not in {os.path.realpath(m) for m in mzmls}:
            os.remove(stale)
            print(f'    removed stale .mzXML symlink: {stale}')

    print(f'    staged {len(staged_paths)} spectra -> {dst}')
    return staged_paths, dst

def copy_annotation_to_stage(annot_path, staged_dir):
    """Copy the (by now fully padded) annotation.txt into an experiment's staged spectra dir --
    this is the file TmtIntegrator actually reads, so it must be the LAST thing written to it,
    strictly after pad_annotation() has finished patching annot_path in annotations/.
    """
    dst = os.path.join(staged_dir, 'annotation.txt')
    shutil.copy(annot_path, dst)
    print(f'    copied padded annotation -> {dst}')

# Files allowed into the final manifest/sample table. .mzXML is deliberately absent: TMTIntegrator
# can't read it, and it should always have been converted to .mzML (or, with --allow-raw, left as
# .raw) well before this point. This is the last line of defense before anything is written out.
ALLOWED_MANIFEST_EXTS = {'.mzml'}

def filter_final_spectra(entries, allow_raw):
    """Keep only .mzML (and .raw if allow_raw) entries; drop + warn about anything else.

    This runs right before write_manifest/write_sample_table so that no .mzXML (or other
    unexpected extension) can ever reach the manifest, regardless of what happened upstream.
    """
    allowed = ALLOWED_MANIFEST_EXTS | ({'.raw'} if allow_raw else set())
    kept, dropped = [], []
    for path, experiment in entries:
        if os.path.splitext(path)[1].lower() in allowed:
            kept.append((path, experiment))
        else:
            dropped.append(path)
    if dropped:
        shown = ', '.join(os.path.basename(p) for p in dropped[:5])
        more = f' ... (+{len(dropped) - 5} more)' if len(dropped) > 5 else ''
        allowed_desc = '.mzML' + (' or .raw' if allow_raw else '')
        print(f'  WARN: dropped {len(dropped)} file(s) from manifest that were not {allowed_desc}: '
              f'{shown}{more}')
    return kept

def patch_line(text, key, value):
    """Replace `key=...` in a workflow, or append it if absent."""
    pat = re.compile(rf'^{re.escape(key)}=.*$', re.MULTILINE)
    line = f'{key}={value}'
    out = pat.sub(line, text)
    return out if out != text else text.rstrip('\n') + f'\n{line}\n'

def write_workflow(template, fasta, channels, workdir, out_path):
    wf = open(template).read()
    wf = patch_line(wf, 'database.db-path', os.path.abspath(fasta))
    wf = patch_line(wf, 'tmtintegrator.channel_num', str(channels))
    # TMT-Integrator writes its report tables to this path -- separate from the general
    # --workdir FragPipe uses for everything else. If it's missing/blank, FragPipe refuses to
    # run with "Output directory can't be left empty". Point it at <workdir>/tmt-report so it's
    # always set and always lands inside this run's own results folder.
    wf = patch_line(wf, 'tmtintegrator.output', os.path.join(os.path.abspath(workdir), 'tmt-report'))
    open(out_path, 'w').write(wf)

def write_manifest(entries, out_path):
    """FragPipe manifest across all experiments: <mzML/raw>\t<experiment>\t1\tDDA."""
    with open(out_path, 'w') as fh:
        for path, experiment in entries:
            fh.write(f'{os.path.abspath(path)}\t{experiment}\t1\tDDA\n')

def write_sample_table(entries, out_path):
    """Human-readable record of file -> experiment (not used by FragPipe)."""
    with open(out_path, 'w') as fh:
        fh.write('file\texperiment\n')
        for path, experiment in entries:
            fh.write(f'{os.path.basename(path)}\t{experiment}\n')

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

def write_submit(run, paths, opts, out_path):
    os.makedirs(paths['logdir'], exist_ok=True)
    os.makedirs(paths['workdir'], exist_ok=True)
    text = SUBMIT_TEMPLATE.format(run=run, partition=opts.partition, threads=opts.threads,
                                  ram=opts.ram, time=opts.time, logdir=paths['logdir'],
                                  java_home=opts.java_home, fragpipe_bin=opts.fragpipe_bin,
                                  workflow=paths['workflow'], manifest=paths['manifest'],
                                  workdir=paths['workdir'], tools_folder=opts.tools_folder)
    open(out_path, 'w').write(text)
    os.chmod(out_path, 0o755)

# Run processing
def main():
    ap = argparse.ArgumentParser(
        description='Prep one FragPipe RUN spanning multiple TMT experiments (plexes) -> one manifest + submit script.')
    ap.add_argument('--run', required=True, help='run token, e.g. cortex_tsumagari (names workflow/manifest/submit + spectra subdir)')
    ap.add_argument('--species', required=True, choices=sorted(SPECIES_TAG))
    ap.add_argument('--experiment', required=True, action='append', nargs=3,
                    metavar=('NAME', 'RAW_DIR', 'PLEX'),
                    help='one TMT experiment: its name, the dir with its .raw/.mzML/.mzXML, and the '
                         'plex key selecting its rows from the sample_map (its TMT plex or sample_ID, '
                         'e.g. 3 or S3). Repeat --experiment for each plex in the run.')
    ap.add_argument('--sample-map', required=True, help='the run\'s multi-plex sample_map .xlsx')
    ap.add_argument('--workflow', required=True, help='FragPipe .workflow template to patch')
    ap.add_argument('--fasta', required=True, help='FragPipe FASTA file to search against (patched into the workflow)')
    ap.add_argument('--out-dir', required=True, help='Frag_outputs root (workflows/manifests/...)')
    ap.add_argument('--spectra-root', required=True)
    ap.add_argument('--trfp', default=os.path.expanduser('~/thermoRawFileParser/ThermoRawFileParser'))
    ap.add_argument('--msconvert', default='msconvert',
                    help="msconvert command for .mzXML -> .mzML (TMTIntegrator cannot read .mzXML "
                         "directly). Quote it if it's more than one word, e.g. "
                         "--msconvert 'singularity exec /scratch/maropakis.a/containers/pwiz.sif msconvert'")
    ap.add_argument('--no-convert', action='store_true', help='skip raw->mzML (spectra already mzML/mzXML)')
    ap.add_argument('--allow-raw', action='store_true',
                    help='keep .raw files as .raw in the manifest instead of converting them to '
                         '.mzML (TMTIntegrator can read .raw natively). Without this flag, .raw is '
                         'always converted to .mzML as before. Either way, .mzXML never reaches the '
                         'manifest -- it is always converted, or dropped with a warning.')
    # submit-script knobs
    ap.add_argument('--fragpipe-bin', default='/home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe')
    ap.add_argument('--tools-folder', default='/home/maropakis.a/fragpipe/fragpipe-24.0/tools')
    ap.add_argument('--java-home', default=os.path.expanduser('~/bin/jdk-17.0.18+8'))
    ap.add_argument('--partition', default='short')
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--ram', type=int, default=64)
    ap.add_argument('--time', default='24:00:00')
    a = ap.parse_args()

    run = a.run.strip().lower()   # token may contain underscores (cortex_tsumagari); keep them
    if not os.path.isfile(a.fasta):
        sys.exit(f'--fasta not found: {a.fasta}')

    wf_dir = os.path.join(a.out_dir, 'workflows')
    mf_dir = os.path.join(a.out_dir, 'manifests')
    annot_dir = os.path.join(a.out_dir, 'annotations')
    submit_dir = os.path.join(a.out_dir, 'submit')
    workdir_path = os.path.join(a.out_dir, 'results', run)  # used for --workdir and tmtintegrator.output
    for d in (wf_dir, mf_dir, annot_dir, submit_dir):
        os.makedirs(d, exist_ok=True)

    print(f'[{run}] species={a.species}  {len(a.experiment)} experiments')

    smap = load_sample_map(a.sample_map)

    # Build every experiment: annotation from its plex rows, then convert + stage its spectra.
    entries = []           # (mzML/raw_path, experiment_name) across all experiments -> single manifest
    channel_counts = {}    # experiment -> count of *labeled* rows in its annotation.txt (diagnostic only)
    plex_sizes = {}        # experiment -> true TMT reagent size detected from filenames (used for channel_num)
    annot_paths = {}       # experiment -> its annotation.txt path (used later to pad to channel_num)
    used_channels = {}     # experiment -> set of physical channels already labeled in its annotation.txt
    staged_dirs = {}       # experiment -> its staged spectra dir (annotation.txt is copied here LAST)
    seen = set()
    for name, raw_dir, plex_key in a.experiment:
        name = name.strip()
        if name in seen:
            sys.exit(f'duplicate experiment name {name!r}')
        seen.add(name)
        if not os.path.isdir(raw_dir):
            sys.exit(f'[{name}] raw_dir not found: {raw_dir}')
        print(f'  experiment {name!r}  plex={plex_key}  raw_dir={raw_dir}')

        sub = select_plex(smap, plex_key)
        annot_path = os.path.join(annot_dir, f'{run}_{name}_annotation.txt')
        channel_counts[name], used_channels[name] = write_annotation(sub, annot_path)
        annot_paths[name] = annot_path
        print(f'    {channel_counts[name]} channels')

        if a.no_convert:
            mzmls = glob_spectra(raw_dir)
            print(f'    --no-convert: using {len(mzmls)} existing .mzML/.mzXML')
        else:
            mzmls = convert_raws(raw_dir, a.trfp, allow_raw=a.allow_raw)
        mzmls = ensure_mzml(mzmls, a.msconvert)  # TMTIntegrator needs .mzML (or .raw), not .mzXML
        if not mzmls:
            sys.exit(f'[{name}]: no mzML/mzXML/raw to stage')

        detected = detect_plex_size(mzmls)
        if detected is not None:
            plex_sizes[name] = detected
        else:
            print(f"    WARN: couldn't detect TMT plex size from filenames in {raw_dir}; "
                  f'will fall back to the {channel_counts[name]} labeled channel(s) for channel_num')

        staged, staged_dir = stage_experiment(mzmls, a.spectra_root, run, name)
        entries.extend((s, name) for s in staged)
        staged_dirs[name] = staged_dir

    # One workflow means one TMT channel_num -- prefer the plex size detected from filenames
    # (the real physical reagent size) over the labeled-channel count, which only reflects how
    # many samples this cohort's sample_map happens to name within a shared/bridged plex.
    if plex_sizes:
        uniq_sizes = set(plex_sizes.values())
        if len(uniq_sizes) > 1:
            print(f'  WARN: experiments were acquired with DIFFERENT TMT plex sizes {plex_sizes} -- '
                  f'one workflow/channel_num cannot correctly quantify all of them together. '
                  f'Split this run by plex size before submitting.')
        channels = max(uniq_sizes)
        undetected = [n for n in channel_counts if n not in plex_sizes]
        if undetected:
            channels = max(channels, max(channel_counts[n] for n in undetected))
    else:
        uniq = set(channel_counts.values())
        if len(uniq) > 1:
            print(f'  WARN: plex size undetectable from filenames and experiments disagree on '
                  f'labeled channel count {channel_counts}; using max (verify this is correct!)')
        channels = max(uniq)

    # channel_num is only known for certain once every experiment has been scanned; now that it
    # is, pad each experiment's annotation.txt with Empty1, Empty2, ... rows for any physical
    # channels the sample_map didn't label, so its line count matches channel_num exactly --
    # THEN (not before) copy that final, padded file into the staged dir TmtIntegrator actually
    # reads from. Doing the copy any earlier is what caused the mismatch error to persist even
    # after pad_annotation() was added: the staged copy would already exist, unpadded, by then.
    for name in seen:
        pad_annotation(annot_paths[name], used_channels[name], channels, name)
        copy_annotation_to_stage(annot_paths[name], staged_dirs[name])

    # Last line of defense: strip anything that isn't .mzML (or .raw, with --allow-raw) before
    # a single line is written to the manifest or sample table.
    entries = filter_final_spectra(entries, a.allow_raw)
    if not entries:
        sys.exit('no .mzML (or .raw) files survived filtering -- nothing to write a manifest for')

    # One FASTA, one workflow, one manifest, one submit for the whole run.
    fasta = a.fasta
    wf_path = os.path.join(wf_dir, f'{run}.workflow')
    mf_path = os.path.join(mf_dir, f'{run}.fp-manifest')
    tbl_path = os.path.join(annot_dir, f'{run}_samples.txt')
    write_workflow(a.workflow, fasta, channels, workdir_path, wf_path)
    write_manifest(entries, mf_path)
    write_sample_table(entries, tbl_path)
    print(f'  channel_num = {channels}\n  workflow -> {wf_path}'
          f'\n  manifest -> {mf_path} ({len(entries)} files, {len(seen)} experiments)'
          f'\n  samples  -> {tbl_path}\n  fasta = {os.path.basename(fasta)}')

    submit_path = os.path.join(submit_dir, f'submit_{run}.sh')
    paths = dict(workflow=wf_path, manifest=mf_path,
                 workdir=workdir_path,
                 logdir=os.path.join(a.out_dir, 'logs'))
    write_submit(run, paths, a, submit_path)
    print(f'  submit   -> {submit_path}\n\nNext: sbatch {submit_path}')


if __name__ == '__main__':
    main()
