#!/usr/bin/env python3 

"""
Script to create per-plex .csv files for buildFragFASTA.py and build correct headers for SAAPs. 

Requires having ran scripts SAAP_Detection & Validation1 from the Decode pipeline in Tsour et al., Nature 2026 

Inputs: 
    *_MTP.fasta - reference FASTA with substituted sequences appended
    
    From MaxQuant Dependent Peptide Search:
        evidence.txt - used for 'Proteins'
        proteinGroups.txt - used for protein accessions, genes, and descriptions 

Outputs:
    {token}.csv - plex-level .csv file schematicized for buildFragFASTA.py:
        sequence (SAAP), accession, gene, description, bp_seq, all_accessions, status, n_base_candidates 

        
Example usage:
  python prepFASTA.py \
    --mtp-dir     /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-root  /scratch/maropakis.a/MQ_outputs/Ping_2018 \
    --human-root  /scratch/maropakis.a/MQ_outputs/Bai_2020 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Takasugi_2024 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Keele_2025 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Tsumagari_2023 \
    --out-dir     /scratch/maropakis.a/Dependencies/mtp_maps/

"""

import argparse
import csv
import os
import re
from collections import defaultdict

# Header parsers 
HDR_RE = re.compile(r'(?:sp|tr)\|([^|]+)\|\S+\s+(.+?)\s+OS=')
GN_RE  = re.compile(r'GN=(\S+)')

# Helper functions
def plex_token(filename): 
    # Function to pull token from *_MTP.fasta file name 
    # e.g. S1_ACGB1_MTP.fasta --> token = 'acgb1'

    stem = re.sub(r'_MTP\.fasta$', '', os.path.basename(filename), flags=re.I)
    stem = re.sub(r'^S\d+_', '', stem)
    return stem.lower()

# Note: depending on data available, some plex tokens get a dataset suffix to avoid tissue-name collisions across studies 
# Add a new disambiguated dataset to DATASET_SUFFIX 
DATASET_SUFFIX = {'keele_2023': 'keele', 'tsumagari_2023': 'tsumagari'}
DATASET_RE = re.compile(r'^([A-Za-z]+(?:_?\d{4}))_(.+)$')   # <Name><year> or <Name>_<year> prefix

def token_from_rest(dataset, rest):
    # Function to build canonical plex token from dataset name + plex/tissue remainder 
    key = dataset.lower()
    if key in DATASET_SUFFIX:
        suf = DATASET_SUFFIX[key]
        if '_' in rest:
            return '_'.join(rest.lower().split('_') + [suf])
        return f'{rest.lower()}_{suf}'
    return re.sub(r'[^A-Za-z0-9]', '', rest).lower()

def dir_token(txt_dir):
    # Function to build token from MQ DP directory names
    # e.g. .../Ping2018_ACG_B1_DP/combined/txt -> 'acgb1' ; Keele2025_kidney_DP -> 'kidney_keele'
    leaf = os.path.basename(os.path.dirname(os.path.dirname(txt_dir)))  # the *_DP folder
    leaf = re.sub(r'_DP$', '', leaf, flags=re.I)
    m = DATASET_RE.match(leaf)
    if not m:
        raise SystemExit(f'cannot parse dataset prefix from DP dir name: {leaf!r}')
    return token_from_rest(m.group(1), m.group(2))

def find_dp_txt_dirs(root):
    # Function to walk a dataset root and return DP combined/txt dirs with evidence.txt + proteinGroups.txt
    txt_dirs = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if 'evidence.txt' in filenames and 'proteinGroups.txt' in filenames:
            leaf = os.path.basename(os.path.dirname(os.path.dirname(dirpath)))
            if leaf.upper().endswith('_DP'):
                txt_dirs.append(dirpath)
    return sorted(txt_dirs)

def build_token_index(roots, species):
    # Function to expand dataset roots into {plex/tissue token: DP txt dir}; error on token collision
    token_to_dir = {}
    for root in roots:
        if not os.path.isdir(root):
            raise SystemExit(f'{species}: root not found: {root}')
        found = find_dp_txt_dirs(root)
        if not found:
            raise SystemExit(f'{species}: no *_DP/combined/txt under {root}')
        print(f'  {species}: {len(found)} DP txt dir(s) under {root}')
        for d in found:
            tok = dir_token(d)
            if tok in token_to_dir:
                raise SystemExit(f'{species}: duplicate DP token {tok!r}:\n'
                                 f'  {token_to_dir[tok]}\n  {d}')
            token_to_dir[tok] = d
    return token_to_dir

def parse_fasta(path):
    # Function to yield (header, sequence); header keeps its leading '>'
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

def parse_protein_groups(mq_dir):
    # Function to parse proteinGroups.txt for accession, gene, and description info
    acc_gene, acc_desc = {}, {}
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

def index_evidence(mq_dir):
    # Function to parse evidence.txt sequences for the 'Proteins' string 
    seq_to_proteins = {}
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
    # Function to extract UniProt accessions in a MaxQuant 'Proteins' string (drops MTP| entries)
    accs = []
    for tok in str(protein_str).split(';'):
        tok = tok.replace('CON__', '')
        if tok.startswith('MTP'):
            continue
        m = re.match(r'(?:sp|tr)\|([^|]+)\|', tok) or re.match(r'^([A-Z0-9][A-Z0-9-]+)$', tok)
        if m:
            accs.append(m.group(1).split('-')[0])
    return list(dict.fromkeys(accs))


# Functions to match SAAP to its BP by a single-residue diference
def build_length_index(sequences):
    # Function to group evidence sequences by length for 0(matches) BP search
    by_len = defaultdict(list)
    for s in sequences:
        by_len[len(s)].append(s)
    return by_len

def find_base_peptides(mtp, by_len):
    # Function to find BP as same-length evidence peptides differing at exactly one position
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

def resolve_one_mtp(mtp, species, by_len, seq_to_proteins, acc_gene, acc_desc):
    # Function to resolve a single SAAP seq against one plex's DP index and return a record dict
    bases = find_base_peptides(mtp, by_len)
    chosen_bp, accs = '', []
    for bp in bases:
        a = proteins_to_accessions(seq_to_proteins.get(bp, ''))
        if a:
            chosen_bp, accs = bp, a
            break
    if accs:
        acc = accs[0]
        return dict(sequence=mtp, species=species, accession=acc,
                    gene=acc_gene.get(acc, ''),
                    description=acc_desc.get(acc, ''),
                    bp_seq=chosen_bp,
                    all_accessions=';'.join(accs),
                    status='keep',
                    n_base_candidates=len(bases))
    return dict(sequence=mtp, species=species, accession='', gene='', description='',
                bp_seq='', all_accessions='', status='unresolved',
                n_base_candidates=len(bases))


# Functions to resolve each SAAP against its own plex/DP dir (important for multiple files in one dir)
def collect_mtp_files(mtp_dir):
    # Function to read all *_MTP.fasta; return [(filename, token, sorted MTP seqs), ...]
    files = []
    seen_tokens = {}
    for fn in sorted(os.listdir(mtp_dir)):
        if not fn.endswith('_MTP.fasta'):
            continue
        tok = plex_token(fn)
        if tok in seen_tokens:
            # why: two MTP files collapsing to one plex token would silently overwrite one CSV.
            raise SystemExit(f'two MTP files share plex token {tok!r}: '
                             f'{seen_tokens[tok]} and {fn}')
        seen_tokens[tok] = fn
        seqs = set()
        for header, seq in parse_fasta(os.path.join(mtp_dir, fn)):
            if header.startswith('>MTP|'):
                s = seq.strip().upper()
                if s:
                    seqs.add(s)
        files.append((fn, tok, sorted(seqs)))
    return files

def write_plex(token, records, out_dir):
    # Function to write one CSV per plex, named '{token}.csv'
    path = os.path.join(out_dir, f'{token}.csv')
    cols = ['sequence', 'species', 'accession', 'gene', 'description',
            'bp_seq', 'all_accessions', 'status', 'n_base_candidates']
    counts = {}
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in records:
            counts[r['status']] = counts.get(r['status'], 0) + 1
            w.writerow(r)
    print(f'  {token}: {len(records)} MTP seqs -> {dict(sorted(counts.items()))}  ({path})')


# Run process
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True,
                    help='dir of *_MTP.fasta (the appended FASTAs)')
    ap.add_argument('--human-root', action='append', default=[],
                    help='dataset root for human (e.g. .../Ping_2018, .../Bai_2020); '
                         'repeatable. Walked for all *_DP/combined/txt dirs.')
    ap.add_argument('--mouse-root', action='append', default=[],
                    help='dataset root for mouse (e.g. .../Takasugi_2024); repeatable. '
                         'Walked for all *_DP/combined/txt dirs.')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    files = collect_mtp_files(a.mtp_dir)
    print(f'Read {len(files)} *_MTP.fasta files ({len(files)} plexes)')

    # one token -> (DP dir, species) index spanning both root sets. species comes from which root
    # flag the dir was found under, so adding a dataset = passing its root, no token rule to edit
    token_index = {}
    for species, roots in (('human', a.human_root), ('mouse', a.mouse_root)):
        if not roots:
            continue
        print(f'{species}:')
        for tok, d in build_token_index(roots, species).items():
            if tok in token_index:
                raise SystemExit(f'plex token {tok!r} found under two roots:\n'
                                 f'  {token_index[tok][0]}\n  {d}')
            token_index[tok] = (d, species)

    # resolve each plex STRICTLY against its own matched DP dir; one CSV per plex, no merging
    print('resolving per plex:')
    for fn, tok, seqs in files:
        if tok not in token_index:
            raise SystemExit(f'{fn}: plex token {tok!r} matched no DP dir. '
                             f'Available tokens: {sorted(token_index)}')
        mq_dir, species = token_index[tok]
        leaf = os.path.basename(os.path.dirname(os.path.dirname(mq_dir)))
        acc_gene, acc_desc = parse_protein_groups(mq_dir)
        seq_to_proteins = index_evidence(mq_dir)
        by_len = build_length_index(seq_to_proteins.keys())
        print(f'  {fn} [{species}] <- {tok} ({leaf}): '
              f'{len(acc_desc):,} proteins, {len(seq_to_proteins):,} peptides, {len(seqs)} MTP seqs')
        records = [resolve_one_mtp(m, species, by_len, seq_to_proteins, acc_gene, acc_desc) for m in seqs]
        write_plex(tok, records, a.out_dir)


if __name__ == '__main__':
    main()