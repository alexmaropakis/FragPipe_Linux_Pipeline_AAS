#!/usr/bin/env python3
"""
7_run_plexes.py 

Per-plex FragPipe workflow + manifest, then headless run.
Consumes 4_stage_spectra.py output: each --spectra-root subfolder is one plex
(lowercase id) holding its spectra + annotation.txt. FASTAs matched
case-insensitively from names like S3_Aorta_fragpipe.fasta.

Example usage:
python3 run_plexes.py --spectra-root /scratch/maropakis.a/spectra \
    --fasta-dir /scratch/maropakis.a/Dependencies/FASTA_fragpipe \
    --template-dir /home/maropakis.a/scripts/FragPipe/templates \
    --out-dir /scratch/maropakis.a/Frag_outputs \
    --fragpipe-bin /home/maropakis.a/fragpipe/fragpipe-24.0/bin/fragpipe \
    --tools-folder /home/maropakis.a/fragpipe/fragpipe-24.0/tools \
    --spectra-ext .mzML

FIXED BUG 06-22-2026: Requires file extension input to prevent cross-talk between .raw and .mzML files in the same folder

"""
import argparse, glob, os, re, subprocess, sys

EXPECTED = {"TMT-10": 10, "TMT-11": 11, "TMT-16": 16}


def route(plex_id):
    tokens = re.split(r"[._\-]", plex_id.upper())
    if any(t.startswith(("FC", "ACG")) for t in tokens):
        return "TMT10_MS3_Val.workflow", "TMT-10"          # human, MS3
    tissues = ("AORTA", "BRAIN", "HEART", "KIDNEY",
               "LIVER", "LUNG", "MUSCLE", "SKIN")
    if any(t.startswith(tissues) for t in tokens):
        return "TMT16_Val.workflow", "TMT-16"              # mouse, MS2
    raise ValueError(f"no template rule (tokens {tokens})")


def patch_line(text, key, value):
    line = f"{key}={value}"
    pat = re.compile(rf"^{re.escape(key)}=.*$", re.MULTILINE)
    out = pat.sub(line, text)
    return out if out != text else text.rstrip("\n") + f"\n{line}\n"


def index_fastas(fasta_dir):
    """lowercased plex label -> path; handles S#_<label>[_MTP]_fragpipe.fasta."""
    idx = {}
    for p in glob.glob(os.path.join(fasta_dir, "*.fasta")):
        label = os.path.basename(p)
        label = re.sub(r"^S\d+_", "", label, flags=re.IGNORECASE)            # drop S#_
        label = re.sub(r"(?:_MTP)?_fragpipe\.fasta$", "", label, flags=re.IGNORECASE)  # drop suffix
        label = re.sub(r"\.fasta$", "", label, flags=re.IGNORECASE)          # any leftover
        idx[label.lower()] = p
    return idx

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra-root", required=True)
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--template-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--fragpipe-bin", required=True)
    ap.add_argument("--tools-folder", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--ram", type=int, default=64)
    ap.add_argument("--only", default=None, help="run just this one plex")
    ap.add_argument("--spectra-ext", default=".mzML",
                    help="spectra file extension to use (e.g. .raw or .mzML)")
    ap.add_argument("--run", action="store_true", help="execute (default: print only)")
    a = ap.parse_args()

    wf_dir = os.path.join(a.out_dir, "workflows")
    mf_dir = os.path.join(a.out_dir, "manifests")
    res_dir = os.path.join(a.out_dir, "results")
    for d in (wf_dir, mf_dir, res_dir):
        os.makedirs(d, exist_ok=True)

    fasta_idx = index_fastas(a.fasta_dir)

    plex_dirs = sorted(d for d in glob.glob(os.path.join(a.spectra_root, "*"))
                       if os.path.isdir(d))
    if not plex_dirs:
        sys.exit(f"No plex folders under {a.spectra_root}")

    for plex_dir in plex_dirs:
        plex = os.path.basename(plex_dir).lower()
        if a.only and plex != a.only.lower():
            continue

        fasta = fasta_idx.get(plex)
        annot = os.path.join(plex_dir, "annotation.txt")
        spectra = sorted(p for p in glob.glob(os.path.join(plex_dir, "*"))
                         if p.lower().endswith(a.spectra_ext.lower()))

        missing = [n for n, ok in [("FASTA", bool(fasta)),
                                   ("annotation.txt", os.path.isfile(annot)),
                                   ("spectra", bool(spectra))] if not ok]
        if missing:
            print(f"SKIP {plex}: missing {', '.join(missing)}")
            continue
        try:
            template, channel = route(plex)
        except ValueError as e:
            print(f"SKIP {plex}: {e}")
            continue

        n_ch = sum(1 for ln in open(annot) if ln.strip())
        if n_ch != EXPECTED[channel]:
            print(f"WARN {plex}: {n_ch} channels but {channel} (expected {EXPECTED[channel]})")

        wf = open(os.path.join(a.template_dir, template)).read()
        wf = patch_line(wf, "database.db-path", os.path.abspath(fasta))
        wf = patch_line(wf, "tmtintegrator.channel_num", channel)
        wf_path = os.path.join(wf_dir, f"{plex}.workflow")
        open(wf_path, "w").write(wf)

        mf_path = os.path.join(mf_dir, f"{plex}.fp-manifest")
        with open(mf_path, "w") as fh:
            for s in spectra:
                fh.write(f"{os.path.abspath(s)}\t{plex}\t1\tDDA\n")

        out = os.path.join(res_dir, plex)
        os.makedirs(out, exist_ok=True)
        cmd = [a.fragpipe_bin, "--headless", "--workflow", wf_path,
               "--manifest", mf_path, "--workdir", out,
               "--threads", str(a.threads), "--ram", str(a.ram),
               "--config-tools-folder", a.tools_folder]
        print(f"\n# {plex}  [{template} / {channel}, {len(spectra)} files, "
              f"{n_ch} channels]  fasta={os.path.basename(fasta)}")
        print(" ".join(cmd))
        if a.run:
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()