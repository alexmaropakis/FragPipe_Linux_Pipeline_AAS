# templates/

![Format](https://img.shields.io/badge/format-.workflow-6E44FF?style=flat-square)
![FragPipe](https://img.shields.io/badge/FragPipe-24.0-6E44FF?style=flat-square)
![Labelling](https://img.shields.io/badge/TMT-10%20%7C%2011%20%7C%2016-F9A825?style=flat-square)
![Quant](https://img.shields.io/badge/quant-MS2%20%7C%20MS3-00838F?style=flat-square)
![Patched](https://img.shields.io/badge/patched-db--path%20%2B%20channel__num-455A64?style=flat-square)

FragPipe `.workflow` templates and an example sample map. The search generators copy one of these
templates and patch only two lines into it — `database.db-path` (the per-plex FASTA) and
`tmtintegrator.channel_num` — so everything else in the file is the fixed search configuration you
choose by picking a template.

Pick a template by TMT plex size, acquisition method, and enzyme.

| Template | Plex | Reporter quant | Enzyme(s) |
|---|---|---|---|
| `TMT10_MS2_Val.workflow` | TMT-10 | MS2 | trypsin |
| `TMT10_MS3_Val.workflow` | TMT-10 | MS3 | trypsin |
| `TMT11_MS2_Val.workflow` | TMT-11 | MS2 | trypsin |
| `TMT11_MS3_Val.workflow` | TMT-11 | MS3 | trypsin |
| `TMT11_MS3_Val_TrypsinLysc.workflow` | TMT-11 | MS3 | trypsin + Lys-C |
| `TMT16_Val.workflow` | TMT-16 | MS2 | trypsin |

Shared across all templates: `rev_` decoy tag, mass calibration on, label-free quant enabled. The
`channel_num` shown is the template default — the generators overwrite it with the actual count
from the sample map, so the plex-size column above is just the intended pairing.

- **MS2 vs MS3** — MS3 (`quant_level=3`) reduces ratio compression at the cost of speed; match it to
  how the data was acquired.
- **Enzyme** — the `_TrypsinLysc` variant adds Lys-C as a second enzyme; the rest are strict trypsin.

## sample_map_example.xlsx

The expected sample-map shape. At minimum the generators need `TMT channel` and `sample_name`
columns (whitespace and case are normalized). For experiment-level runs, add a plex column
(`TMT plex` or `sample_ID`) so each plex's rows can be selected by key.
