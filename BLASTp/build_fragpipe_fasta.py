#!/usr/bin/env python3
"""
Build FragPipe-ready FASTAs from *_MTP.fasta using the per-species filtered CSVs.
Reference entries pass through unchanged; kept MTP entries get a Philosopher-safe
mock-UniProt header whose accession is acc-<mtp_id>-<tmt_set> (unique per entry);
reversed rev_ decoys are appended for every target.
ACG*/FC* -> human; everything else -> mouse.

  python3 build_fragpipe_fasta.py \
    --mtp-dir   /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-csv /scratch/maropakis.a/Dependencies/mtp_maps/human_filtered.csv \
    --mouse-csv /scratch/maropakis.a/Dependencies/mtp_maps/mouse_filtered.csv \
    --out-dir   /scratch/maropakis.a/Dependencies/FASTA_fragpipe/
"""
import argparse, csv, os, re

SPECIES_TAG = {'human': ('Homo sapiens', 9606),
               'mouse': ('Mus musculus', 10090)}
MTP_ID_RE = re.compile(r'(MTP\d+)\b')

def species_for(filename):
    tokens = re.split(r'[._\-]', os.path.basename(filename).upper())
    is_human = any(t.startswith('ACG') or t.startswith('FC') for t in tokens)
    return 'human' if is_human else 'mouse'

def tmt_set_for(filename):
    # MUST match blast_mtps.py: 'ACG_B5_MTP.fasta' -> 'acgb5'
    stem = re.sub(r'_MTP\.fasta$', '', os.path.basename(filename), flags=re.I)
    return re.sub(r'[^A-Za-z0-9]', '', stem).lower()

def load_keep(csv_path):
    """sequence -> (accession, gene) for kept MTPs (classification is seq-determined)."""
    keep = {}
    if csv_path and os.path.exists(csv_path):
        for row in csv.DictReader(open(csv_path)):
            if row['status'] == 'keep':
                keep[row['sequence']] = (row['accession'], row['gene'])
    return keep

def parse_fasta(path):
    ent, h, s = [], None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith('>'):
            if h is not None: ent.append((h, ''.join(s)))
            h, s = line, []
        else: s.append(line)
    if h is not None: ent.append((h, ''.join(s)))
    return ent

def mtp_header(accession, gene, mid, tmt, species):
    os_name, ox = SPECIES_TAG[species]
    return (f'>sp|{accession}|{gene}-mut {gene} mistranslated {mid} '
            f'TMT={tmt} OS={os_name} OX={ox} GN={gene} PE=1 SV=1')

def build(src, dst, keep, species, prefix='rev_'):
    tmt = tmt_set_for(src)
    targets, seen_acc, n_ref, n_mtp = [], set(), 0, 0
    for h, seq in parse_fasta(src):
        if h.startswith('>MTP|'):
            if seq in keep:
                acc, gene = keep[seq]
                m = MTP_ID_RE.search(h)
                mid = m.group(1) if m else h[1:].split()[0]
                accession, base, k = f'{acc}-{mid}-{tmt}', f'{acc}-{mid}-{tmt}', 2
                while accession in seen_acc:        # guard residual dup within file
                    accession = f'{base}-d{k}'; k += 1
                seen_acc.add(accession)
                targets.append((mtp_header(accession, gene, mid, tmt, species), seq))
                n_mtp += 1
        else:
            targets.append((h, seq)); n_ref += 1     # reference, untouched
    with open(dst, 'w') as out:
        for h, seq in targets:
            out.write(h + '\n' + seq + '\n')
        for h, seq in targets:
            out.write('>' + prefix + h[1:] + '\n' + seq[::-1] + '\n')
    return n_ref, n_mtp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True)
    ap.add_argument('--human-csv', required=True)
    ap.add_argument('--mouse-csv', required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    keep = {'human': load_keep(a.human_csv), 'mouse': load_keep(a.mouse_csv)}

    n = 0
    for fn in sorted(os.listdir(a.mtp_dir)):
        if not fn.endswith('_MTP.fasta'): continue
        sp = species_for(fn)
        dst = os.path.join(a.out_dir, fn.replace('_MTP.fasta', '_fragpipe.fasta'))
        n_ref, n_mtp = build(os.path.join(a.mtp_dir, fn), dst, keep[sp], sp)
        print(f'{fn} [{sp}]: {n_ref} ref + {n_mtp} MTP kept')
        n += 1
    print(f'Wrote {n} FASTAs')

if __name__ == '__main__':
    main()
