#!/usr/bin/env python3
"""
1_PrepFASTA.py 

Create .csv files for 2_build_FragFASTA.py to use and build correct headers for MTPs 

SAAP list is the appended *_MTP.fasta, parent proteins are derived from MaxQuant DP search 
txt outputs:
    BP_seq - evidence.txt --> 'Proteins' 
    --> accession - proteinGroups.txt 'Fasta headers' (GN=) 
    --> gene, description 

*_MTP.fasta only has the peptide sequence, but BP is normal MQ identification, so it appears 
in the 'Sequence' row; each MTP is matched to its BP by a single-residue difference (same length,
one mismatch) within the indexed evidence seqs, then resolved as above.

Every MTP whose base peptide resolves to a real accession 
is kept (status='keep'); the rest are 'unresolved'.

Emits {species}.csv with the schema 2_build_FragFASTA.py consumes:
    sequence , accession , gene , description , bp_seq , all_accessions , status
  - `sequence` is the MTP peptide (build_fragpipe_fasta.py keys keep[seq] on it)

Species routing mirrors the pipeline: file tokens starting with ACG/FC -> human;
everything else -> mouse.

Example usage: 
-------
  python3 1_PrepFASTA.py \
    --mtp-dir     /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-root  /scratch/maropakis.a/MQ_outputs/Ping_2018 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Takasugi_2024 \
    --out-dir     /scratch/maropakis.a/Dependencies/mtp_maps/

Each root is walked recursively for every <plex_or_tissue>/combined/txt/ dir
holding evidence.txt + proteinGroups.txt (e.g. Ping_2018/acg/b1/combined/txt,
Ping_2018/fc/.../combined/txt, Takasugi_2024/aorta/combined/txt). Species is
decided by which root flag the dataset is passed under, not by path tokens.
"""

import argparse
import csv
import os
import re
from collections import defaultdict

## header parsers 
HDR_RE = re.compile(r'(?:sp|tr)\|([^|]+)\|\S+\s+(.+?)\s+OS=')
GN_RE  = re.compile(r'GN=(\S+)')

## Helper functions 
def species_for(name):
    """Route a *_MTP.fasta filename to a species"""
    tokens = re.split(r'[._\-\s]', os.path.basename(str(name)).upper())
    is_human = any(t.startswith('ACG') or t.startswith('FC') for t in tokens)
    return 'human' if is_human else 'mouse'

def find_txt_dirs(root):
    """Walk a dataset root; return every dir holding evidence.txt + proteinGroups.txt.

    MQ outputs are nested as <root>/<plex_or_tissue>/combined/txt/, possibly more
    than one level deep (e.g. Ping_2018/acg/b1/combined/txt). This finds them all
    regardless of intermediate structure.
    """
    txt_dirs = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if 'evidence.txt' in filenames and 'proteinGroups.txt' in filenames:
            txt_dirs.append(dirpath)
    return sorted(txt_dirs)

def collect_txt_dirs(roots, species):
    """Expand a list of dataset roots into all their combined/txt dirs."""
    txt_dirs = []
    for root in roots:
        if not os.path.isdir(root):
            raise SystemExit(f'{species}: root not found: {root}')
        found = find_txt_dirs(root)
        if not found:
            raise SystemExit(f'{species}: no evidence.txt+proteinGroups.txt under {root}')
        print(f'  {species}: {len(found)} MQ txt dir(s) under {root}')
        txt_dirs.extend(found)
    return txt_dirs

def parse_fasta(path):
    """Yield (header, sequence); header keeps its leading '>'."""
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


## proteinGroups.txt: accession -> (gene, description)
def parse_protein_groups(mq_dirs):
    acc_gene, acc_desc = {}, {}
    for mq_dir in mq_dirs:
        pg = os.path.join(mq_dir, 'proteinGroups.txt')
        if not os.path.exists(pg):
            raise SystemExit(f'proteinGroups.txt not found: {pg}')
        with open(pg, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            if 'Fasta headers' not in (reader.fieldnames or []):
                raise SystemExit(f"'Fasta headers' missing in {pg}. "
                                 f'Available: {reader.fieldnames}')
            for row in reader:
                for entry in (row.get('Fasta headers') or '').split(';'):
                    m = HDR_RE.search(entry)
                    if not m:
                        continue
                    acc_desc.setdefault(m.group(1), m.group(2).strip())
                    g = GN_RE.search(entry)
                    if g:
                        acc_gene.setdefault(m.group(1), g.group(1))
    return acc_gene, acc_desc


## evidence.txt: Sequence -> 'Proteins' string 
def index_evidence(mq_dirs):
    seq_to_proteins = {}
    for mq_dir in mq_dirs:
        ev = os.path.join(mq_dir, 'evidence.txt')
        if not os.path.exists(ev):
            raise SystemExit(f'evidence.txt not found: {ev}')
        with open(ev, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            cols = reader.fieldnames or []
            if 'Sequence' not in cols or 'Proteins' not in cols:
                raise SystemExit(f"'Sequence'/'Proteins' missing in {ev}. "
                                 f'Available: {cols}')
            for row in reader:
                seq = (row.get('Sequence') or '').strip().upper()
                prot = (row.get('Proteins') or '').strip()
                if seq and prot and prot.lower() != 'nan':
                    seq_to_proteins.setdefault(seq, prot)
    return seq_to_proteins

def proteins_to_accessions(protein_str):
    """Real UniProt accessions in a MaxQuant 'Proteins' string (drops MTP| entries)."""
    accs = []
    for tok in str(protein_str).split(';'):
        tok = tok.replace('CON__', '')
        if tok.startswith('MTP'):
            continue
        m = re.match(r'(?:sp|tr)\|([^|]+)\|', tok) or re.match(r'^([A-Z0-9][A-Z0-9-]+)$', tok)
        if m:
            accs.append(m.group(1).split('-')[0])
    return list(dict.fromkeys(accs))


## MTP -> base peptide by single-residue difference 
def build_length_index(sequences):
    """Group evidence sequences by length for O(matches) base-peptide search."""
    by_len = defaultdict(list)
    for s in sequences:
        by_len[len(s)].append(s)
    return by_len

def find_base_peptides(mtp, by_len):
    """Base candidates = same-length evidence peptides differing at exactly one position."""
    out = []
    for cand in by_len.get(len(mtp), ()):
        if cand == mtp:
            continue
        diff = 0
        for a, b in zip(mtp, cand):
            if a != b:
                diff += 1
                if diff > 1:
                    break
        if diff == 1:
            out.append(cand)
    return out

def resolve_species(mtp_seqs, mq_dirs):
    acc_gene, acc_desc = parse_protein_groups(mq_dirs)
    print(f'  {len(acc_desc):,} reference proteins parsed from FASTA headers')
    seq_to_proteins = index_evidence(mq_dirs)
    print(f'  {len(seq_to_proteins):,} unique peptide sequences indexed '
          f'across {len(mq_dirs)} evidence.txt file(s)')
    by_len = build_length_index(seq_to_proteins.keys())

    records = []
    for mtp in mtp_seqs:
        bases = find_base_peptides(mtp, by_len)
        # pick the first base whose Proteins resolves to a real accession
        chosen_bp, accs = '', []
        for bp in bases:
            a = proteins_to_accessions(seq_to_proteins.get(bp, ''))
            if a:
                chosen_bp, accs = bp, a
                break
        if accs:
            acc = accs[0]
            records.append(dict(sequence=mtp, accession=acc,
                                gene=acc_gene.get(acc, ''),
                                description=acc_desc.get(acc, ''),
                                bp_seq=chosen_bp,
                                all_accessions=';'.join(accs),
                                status='keep',
                                n_base_candidates=len(bases)))
        else:
            records.append(dict(sequence=mtp, accession='', gene='',
                                description='', bp_seq='', all_accessions='',
                                status='unresolved', n_base_candidates=len(bases)))
    return records

def write_species(species, records, out_dir):
    path = os.path.join(out_dir, f'{species}.csv')
    cols = ['sequence', 'accession', 'gene', 'description',
            'bp_seq', 'all_accessions', 'status', 'n_base_candidates']
    counts = {}
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in records:
            counts[r['status']] = counts.get(r['status'], 0) + 1
            w.writerow(r)
    print(f'{species}: {len(records)} MTP seqs -> {dict(sorted(counts.items()))}  ({path})')

def collect_mtp_seqs(mtp_dir):
    """Read all *_MTP.fasta; return {species: set(MTP peptide sequences)}."""
    groups = {'human': set(), 'mouse': set()}
    n_files = 0
    for fn in sorted(os.listdir(mtp_dir)):
        if not fn.endswith('_MTP.fasta'):
            continue
        sp = species_for(fn)
        for header, seq in parse_fasta(os.path.join(mtp_dir, fn)):
            if header.startswith('>MTP|'):
                s = seq.strip().upper()
                if s:
                    groups[sp].add(s)
        n_files += 1
    return groups, n_files

## Run 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True,
                    help='dir of *_MTP.fasta (the appended FASTAs)')
    ap.add_argument('--human-root', action='append', default=[],
                    help='dataset root for human (e.g. .../Ping_2018); repeatable. '
                         'Walked for all combined/txt dirs.')
    ap.add_argument('--mouse-root', action='append', default=[],
                    help='dataset root for mouse (e.g. .../Takasugi_2024); repeatable. '
                         'Walked for all combined/txt dirs.')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    groups, n_files = collect_mtp_seqs(a.mtp_dir)
    print(f'Read {n_files} *_MTP.fasta files: '
          f"{len(groups['human'])} human, {len(groups['mouse'])} mouse unique MTP seqs")

    for species, roots in (('human', a.human_root), ('mouse', a.mouse_root)):
        seqs = sorted(groups[species])
        if not seqs:
            print(f'{species}: no MTP seqs, skipping')
            continue
        if not roots:
            raise SystemExit(f'{species}: {len(seqs)} MTP seqs but no --{species}-root given')
        print(f'{species}:')
        mq_dirs = collect_txt_dirs(roots, species)
        records = resolve_species(seqs, mq_dirs)
        write_species(species, records, a.out_dir)


if __name__ == '__main__':
    main()