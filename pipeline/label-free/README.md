# label-free/

![Status](https://img.shields.io/badge/status-placeholder-9E9E9E?style=flat-square)
![Quant](https://img.shields.io/badge/quant-label--free-lightgrey?style=flat-square)
![Implemented](https://img.shields.io/badge/implemented-not%20yet-D84315?style=flat-square)

Placeholder for a label-free (LFQ) variant of the pipeline. Nothing here yet.

The `per-plex/` and `experiment-level-plex/` generators are TMT-only — channel order, annotation
files, and `tmtintegrator.channel_num` are all TMT-specific. A label-free generator would drop the
reporter-channel machinery: no annotation.txt, no `channel_num`, and IonQuant/FreeQuant LFQ in the
workflow instead of TMT-Integrator. The FASTA build (`utility/`) and raw → mzML conversion carry
over unchanged.
