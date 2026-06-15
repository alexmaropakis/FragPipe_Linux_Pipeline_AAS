#!/usr/bin/env python3
"""
2_buildFragFASTA.py 

Build FragPipe-ready FASTAs from *_MTP.fasta using the per-species filtered CSVs.

Reference entries pass through unchanged. Kept MTP entries get a Philosopher-safe
mock-UniProt header whose accession is <acc>-<mtp_id>-<tmt_set> (made unique per
entry). Reversed rev_ decoys are appended for *every* target, including MTPs.

ACG*/FC* -> human; everything else -> mouse.

Example:
  python3 2_buildFragFASTA.py \
    --mtp-dir   /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-csv /scratch/maropakis.a/Dependencies/mtp_maps/human.csv \
    --mouse-csv /scratch/maropakis.a/Dependencies/mtp_maps/mouse.csv \
    --out-dir   /scratch/maropakis.a/Dependencies/FASTA_fragpipe/

Run after 1_PrepFASTA.py

"""

import argparse
import csv
import os
import re

SPECIES_TAG = {
    'human': ('Homo sapiens', 9606),
    'mouse': ('Mus musculus', 10090),
}
MTP_ID_RE = re.compile(r'(MTP\d+)\b')
DECOY_PREFIX = 'rev_'

## Helper functions 
def species_for(filename):
    tokens = re.split(r'[._\-]', os.path.basename(filename).upper())
    is_human = any(t.startswith('ACG') or t.startswith('FC') for t in tokens)
    return 'human' if is_human else 'mouse'


def tmt_set_for(filename):
    """e.g. 'ACG_B5_MTP.fasta' -> 'acgb5'."""
    stem = re.sub(r'_MTP\.fasta$', '', os.path.basename(filename), flags=re.I)
    return re.sub(r'[^A-Za-z0-9]', '', stem).lower()


def load_keep(csv_path):
    """Return {sequence: (accession, gene)} for MTPs marked 'keep'."""
    keep = {}
    if csv_path and os.path.exists(csv_path):
        with open(csv_path, newline='') as fh:
            for row in csv.DictReader(fh):
                if row['status'] == 'keep':
                    keep[row['sequence']] = (row['accession'], row['gene'])
    return keep

def parse_fasta(path):
    """Yield (header, sequence) tuples; header keeps its leading '>'."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header, seq = line, []
            else:
                seq.append(line)
    if header is not None:
        yield header, ''.join(seq)

def mtp_header(accession, gene, mid, species):
    os_name, ox = SPECIES_TAG[species]
    return (f'>sp|{accession}|{gene}-mut {gene} mistranslated {mid} '
            f'OS={os_name} OX={ox} GN={gene} PE=1 SV=1')

def unique_accession(base, seen):
    """Return base, or base-d2, base-d3, ... if already taken."""
    accession, k = base, 2
    while accession in seen:
        accession = f'{base}-d{k}'
        k += 1
    seen.add(accession)
    return accession

def build(src, dst, keep, species):
    tmt = tmt_set_for(src)
    targets, seen_acc = [], set()
    n_ref = n_mtp = 0

    for header, seq in parse_fasta(src):
        if header.startswith('>MTP|'):
            if seq not in keep:
                continue
            acc, gene = keep[seq]
            m = MTP_ID_RE.search(header)
            mid = m.group(1) if m else header[1:].split()[0]
            accession = unique_accession(f'{acc}-{mid}-{tmt}', seen_acc)
            targets.append((mtp_header(accession, gene, mid, species), seq))
            n_mtp += 1
        else:
            targets.append((header, seq))
            n_ref += 1

    with open(dst, 'w') as out:
        # Forward targets.
        for header, seq in targets:
            out.write(f'{header}\n{seq}\n')
        # Reversed decoys for every target (references AND MTPs).
        for header, seq in targets:
            out.write(f'>{DECOY_PREFIX}{header[1:]}\n{seq[::-1]}\n')

    return n_ref, n_mtp

## Run 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True)
    ap.add_argument('--human-csv', required=True)
    ap.add_argument('--mouse-csv', required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    keep = {
        'human': load_keep(a.human_csv),
        'mouse': load_keep(a.mouse_csv),
    }

    n = 0
    for fn in sorted(os.listdir(a.mtp_dir)):
        if not fn.endswith('_MTP.fasta'):
            continue
        sp = species_for(fn)
        dst = os.path.join(a.out_dir, fn.replace('_MTP.fasta', '_fragpipe.fasta'))
        n_ref, n_mtp = build(os.path.join(a.mtp_dir, fn), dst, keep[sp], sp)
        print(f'{fn} [{sp}]: {n_ref} ref + {n_mtp} MTP kept')
        n += 1

    print(f'Wrote {n} FASTAs')


if __name__ == '__main__':
    main()
