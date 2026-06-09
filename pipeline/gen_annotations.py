#!/usr/bin/env python3
# make_annotations.py — one FragPipe annotation.txt per TMT set from sample_map xlsx
import glob, os, re, pandas as pd
from collections import Counter

MAP_DIR = '/scratch/maropakis.a/Dependencies/sample_map/'
OUT_DIR = '/scratch/maropakis.a/Dependencies/annotations/'
os.makedirs(OUT_DIR, exist_ok=True)

# canonical channel order; present channels are sorted into this, gaps are fine
ORDER = ['126','127N','127C','128N','128C','129N','129C','130N','130C',
         '131','131N','131C','132N','132C','133N','133C','134N']
ORD = {c: i for i, c in enumerate(ORDER)}

def norm(col):                       # 'Sample name'/'sample_name'/'TMT channel' -> snake
    return re.sub(r'\s+', '_', str(col).strip().lower())

def set_slug(path):                  # sample_map_acgb5.xlsx -> acgb5
    stem = os.path.splitext(os.path.basename(path))[0]
    return re.sub(r'^sample_map_', '', stem, flags=re.I).lower()

written = 0
for path in sorted(glob.glob(os.path.join(MAP_DIR, '*.xlsx'))):
    df = pd.read_excel(path)
    df.columns = [norm(c) for c in df.columns]
    if not {'tmt_channel', 'sample_name'} <= set(df.columns):
        print(f'SKIP {os.path.basename(path)}: no tmt_channel/sample_name column')
        continue

    df = df.dropna(subset=['tmt_channel', 'sample_name'])
    df['tmt_channel'] = df['tmt_channel'].astype(str).str.strip()
    df['sample_name'] = df['sample_name'].astype(str).str.strip()

    # disambiguate labels that repeat within the plex (e.g. GIS on 126 and 131)
    counts = Counter(df['sample_name'])
    rows = [(ch, f'{name}_{ch}' if counts[name] > 1 else name)
            for ch, name in zip(df['tmt_channel'], df['sample_name'])]
    rows.sort(key=lambda x: ORD.get(x[0], 999))

    slug = set_slug(path)
    out = os.path.join(OUT_DIR, f'{slug}_annotation.txt')
    with open(out, 'w') as f:
        for ch, name in rows:
            f.write(f'{ch} {name}\n')
    written += 1
    print(f'{slug}: {len(rows)} ch -> {out}')

print(f'\nWrote {written} annotation files to {OUT_DIR}')
