#!/usr/bin/env python3
"""
Map MTP/SAAP peptides to parent proteins via blastp, per species.
ACG*/FC* *_MTP.fasta -> human; everything else -> mouse.
Writes 4 CSVs: {species}_unfiltered.csv and {species}_filtered.csv.
Requires NCBI BLAST+ on PATH.

Query ids are the FULL MTP header + TMT-set tag (e.g. '...|TMT=acgb5'), so the
same MTP number reused across TMT sets never overwrites another in the query
FASTA or in the per-query hit grouping.

  python3 blast_mtps.py \
    --mtp-dir    /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-ref  /scratch/maropakis.a/Dependencies/FASTA/HUMAN.fasta \
    --mouse-ref  /scratch/maropakis.a/Dependencies/FASTA/MOUSE_UP000000589_10090.fasta \
    --out-dir    /scratch/maropakis.a/Dependencies/mtp_maps/ \
    --threads 16
"""
import argparse, csv, os, re, subprocess

IG_GENE    = re.compile(r'^(IG[HKL]|JCHAIN|IGJ)', re.I)
TRYP_GENE  = re.compile(r'^(PRSS[123]|TRY\d*)$', re.I)
IG_TITLE   = re.compile(r'immunoglobulin', re.I)
TRYP_TITLE = re.compile(r'trypsin', re.I)
GN_RE      = re.compile(r'\bGN=(\S+)')
MTP_ID_RE  = re.compile(r'(MTP\d+)\b')

def species_for(filename):
    # Human samples are tokens beginning with ACG or FC (e.g. ACGB1, FCB3);
    # everything else (Aorta, Brain, Heart, ...) is mouse.
    tokens = re.split(r'[._\-]', os.path.basename(filename).upper())
    is_human = any(t.startswith('ACG') or t.startswith('FC') for t in tokens)
    return 'human' if is_human else 'mouse'

def tmt_set_for(filename):
    # TMT-set tag from the filename, e.g. 'ACG_B5_MTP.fasta' -> 'acgb5'.
    # Tweak this if you want a different slug (keep it whitespace-free).
    stem = re.sub(r'_MTP\.fasta$', '', os.path.basename(filename), flags=re.I)
    return re.sub(r'[^A-Za-z0-9]', '', stem).lower()

def parse_mtp_entries(path):
    """Yield (full_header_without_'>', sequence) for every >MTP| record."""
    h, s = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith('>'):
            if h and h.startswith('>MTP|'):
                yield h[1:], ''.join(s)
            h, s = line, []
        else:
            s.append(line)
    if h and h.startswith('>MTP|'):
        yield h[1:], ''.join(s)

def collect_records(files):
    """One record per MTP entry; query id = full header + TMT tag, made unique."""
    records, seen = [], {}                  # seen: qid -> seq (dup guard)
    for path in files:
        tmt = tmt_set_for(path)
        for raw_header, seq in parse_mtp_entries(path):
            if not seq:
                continue
            token = re.sub(r'\s+', '_', raw_header)        # whole header, no whitespace
            qid   = f'{token}|TMT={tmt}'
            if qid in seen:
                if seen[qid] == seq:
                    continue                                # exact duplicate -> skip
                k = 2
                while f'{qid}|dup{k}' in seen:
                    k += 1
                qid = f'{qid}|dup{k}'                        # same id, diff seq -> uniquify
            seen[qid] = seq
            m = MTP_ID_RE.search(raw_header)
            mtp_id = m.group(1) if m else raw_header.split()[0]
            records.append(dict(qid=qid, mtp_id=mtp_id, tmt=tmt, seq=seq,
                                source=os.path.basename(path)))
    return records

def ensure_db(ref):
    if not os.path.exists(ref + '.phr'):
        subprocess.run(['makeblastdb', '-in', ref, '-dbtype', 'prot',
                        '-parse_seqids'], check=True)

def run_blast(query, ref, out_tab, threads):
    fmt = ('6 qseqid sseqid pident length mismatch gapopen '
           'qstart qend qlen stitle bitscore')
    subprocess.run(['blastp', '-task', 'blastp-short', '-query', query,
                    '-db', ref, '-outfmt', fmt, '-out', out_tab,
                    '-evalue', '1000', '-comp_based_stats', '0',
                    '-max_target_seqs', '50', '-num_threads', str(threads)],
                   check=True)

def acc_gene(sseqid, stitle):
    acc, entry = sseqid, ''
    parts = sseqid.split('|')
    if len(parts) >= 3:                 # sp|P04637|P53_HUMAN
        acc, entry = parts[1], parts[2]
    m = GN_RE.search(stitle)
    gene = m.group(1) if m else (entry.split('_')[0] if entry else acc)
    return acc, gene

def classify(hits):
    full = [h for h in hits if h['gapopen'] == 0 and h['qstart'] == 1
            and h['qend'] == h['qlen'] and h['length'] == h['qlen']]
    if any(h['mismatch'] == 0 for h in full):
        return {'status': 'drop_reference'}
    cand = sorted((h for h in full if h['mismatch'] == 1),
                  key=lambda h: -h['bitscore'])
    if not cand:
        return {'status': 'drop_no_parent'}
    best = cand[0]
    a, g = acc_gene(best['sseqid'], best['stitle'])
    if TRYP_GENE.match(g) or TRYP_TITLE.search(best['stitle']):
        return {'status': 'drop_trypsin', 'accession': a, 'gene': g}
    if IG_GENE.match(g) or IG_TITLE.search(best['stitle']):
        return {'status': 'drop_ig', 'accession': a, 'gene': g}
    genes = {acc_gene(h['sseqid'], h['stitle'])[1] for h in cand}
    return {'status': 'keep', 'accession': a, 'gene': g,
            'ambiguous': len(genes) > 1}

def process_species(species, files, ref, out_dir, threads):
    if not files:
        print(f'{species}: no files, skipping'); return
    records = collect_records(files)
    query = os.path.join(out_dir, f'mtp_query_{species}.fasta')
    tab   = os.path.join(out_dir, f'mtp_blast_{species}.tsv')
    with open(query, 'w') as f:
        for r in records:
            f.write(f">{r['qid']}\n{r['seq']}\n")
    ensure_db(ref)
    run_blast(query, ref, tab, threads)

    by_q = {}
    for line in open(tab):
        p = line.rstrip('\n').split('\t')
        by_q.setdefault(p[0], []).append(dict(
            sseqid=p[1], length=int(p[3]), mismatch=int(p[4]),
            gapopen=int(p[5]), qstart=int(p[6]), qend=int(p[7]),
            qlen=int(p[8]), stitle=p[9], bitscore=float(p[10])))

    unfilt = os.path.join(out_dir, f'{species}_unfiltered.csv')
    filt   = os.path.join(out_dir, f'{species}_filtered.csv')
    counts = {}
    cols = ['query_id', 'mtp_id', 'tmt_set', 'source_file', 'sequence',
            'species', 'status', 'accession', 'gene', 'ambiguous']
    with open(unfilt, 'w', newline='') as uf, open(filt, 'w', newline='') as ff:
        wu, wf = csv.writer(uf), csv.writer(ff)
        wu.writerow(cols); wf.writerow(cols)
        for r in records:
            res = classify(by_q.get(r['qid'], []))
            counts[res['status']] = counts.get(res['status'], 0) + 1
            row = [r['qid'], r['mtp_id'], r['tmt'], r['source'], r['seq'],
                   species, res['status'], res.get('accession', ''),
                   res.get('gene', ''), res.get('ambiguous', '')]
            wu.writerow(row)
            if res['status'] == 'keep':
                wf.writerow(row)
    uniq = len({r['seq'] for r in records})
    print(f'{species}: {len(records)} MTP entries ({uniq} unique seqs) ->',
          dict(sorted(counts.items())))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtp-dir', required=True)
    ap.add_argument('--human-ref', required=True)
    ap.add_argument('--mouse-ref', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    groups = {'human': [], 'mouse': []}
    for fn in sorted(os.listdir(a.mtp_dir)):
        if fn.endswith('_MTP.fasta'):
            groups[species_for(fn)].append(os.path.join(a.mtp_dir, fn))

    process_species('human', groups['human'], a.human_ref, a.out_dir, a.threads)
    process_species('mouse', groups['mouse'], a.mouse_ref, a.out_dir, a.threads)

if __name__ == '__main__':
    main()
