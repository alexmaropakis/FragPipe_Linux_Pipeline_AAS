#!/usr/bin/env python3
"""stage_spectra.py — build one-folder-per-plex (symlinks) + annotation."""
import argparse, glob, os, shutil

EXTS = (".raw", ".mzML", ".mzml")
ONE_PLEX_PER_FOLDER = {
    "Takasugi_2024": ("aorta", "brain", "heart", "kidney",
                      "liver", "lung", "muscle", "skin"),
}
NESTED_BATCH = {"Ping_2018": {"ACG": "acg", "FC": "fc"}}  # group/bN -> prefix+bN


def stage(files, plex, annot_dir, spectra_root):
    dst = os.path.join(spectra_root, plex)
    os.makedirs(dst, exist_ok=True)
    for f in files:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
    ann = os.path.join(annot_dir, f"{plex}_annotation.txt")
    if os.path.isfile(ann):
        shutil.copy(ann, os.path.join(dst, "annotation.txt"))
        print(f"  staged {plex}: {len(files)} files + annotation")
    else:
        print(f"  WARN {plex}: staged {len(files)} files but NO annotation ({ann})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-root", required=True)
    ap.add_argument("--annot-dir", required=True)
    ap.add_argument("--spectra-root", required=True)
    a = ap.parse_args()
    os.makedirs(a.spectra_root, exist_ok=True)

    for dataset, plexes in ONE_PLEX_PER_FOLDER.items():
        for plex in plexes:
            folder = os.path.join(a.raw_root, dataset, plex)
            if not os.path.isdir(folder):
                continue
            files = [p for p in glob.glob(os.path.join(folder, "*")) if p.endswith(EXTS)]
            stage(files, plex, a.annot_dir, a.spectra_root)

    for dataset, groups in NESTED_BATCH.items():
        for group, prefix in groups.items():
            group_dir = os.path.join(a.raw_root, dataset, group)
            if not os.path.isdir(group_dir):
                continue
            for batch in sorted(os.listdir(group_dir)):
                batch_dir = os.path.join(group_dir, batch)
                if not os.path.isdir(batch_dir):
                    continue
                files = [p for p in glob.glob(os.path.join(batch_dir, "*")) if p.endswith(EXTS)]
                if files:
                    stage(files, f"{prefix}{batch}", a.annot_dir, a.spectra_root)


if __name__ == "__main__":
    main()

""" 
python3 stage_spectra.py \
  --raw-root     /scratch/maropakis.a/MQ_raw \
  --annot-dir    /scratch/maropakis.a/Dependencies/annotations \
  --spectra-root /scratch/maropakis.a/spectra
"""
