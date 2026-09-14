# figures/ — authoritative figure set

These files follow the numbering of the manuscript (BMC Genomics submission):
**Fig1–Fig8** are the main-text figures and **FigS1** belongs to Additional file 1
(the exploratory diagnostic scheme). Every figure is provided as PNG (300 dpi) and
PDF; figure captions are the accompanying `*_caption.txt` files.

The files here are byte-identical to the finalized figure set of the submission.

## What changed, and why

An earlier generation of this repository carried `Figure1`–`Figure8` under a scheme in
which the diagnostic-scheme figure sat at position 2, shifting every later number by one.
Those files, and their contents, are superseded; they have been removed. `Figure8.svg`
was the only vector original of the old Fig. 8 and is also gone — the PDFs here are vector.

## How to regenerate

`figure_scripts/` holds the five data-driven scripts that produced Fig. 4–8; every value
they plot is computed from the source tables at run time (no hard-coded numbers). See
`figure_scripts/README_数据来源.md` for the script-to-data mapping.

> ⚠️ `supplementary/FigureS1`–`FigureS8` elsewhere in this repository are an **older**
> supplementary-figure scheme and do **not** correspond to the manuscript's current
> Additional file 1: Fig. S1–S5. Treat `figures/` as authoritative for the main text.
