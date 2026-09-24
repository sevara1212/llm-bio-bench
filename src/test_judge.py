"""Sanity check for judge.py: tricky answer strings with a known correct label.

Usage: python src/test_judge.py   (~20 API calls, a few cents)
Rerun whenever INSTRUCTIONS, DESCRIBE or the label lists in judge.py change.
"""
from judge import NOT_PBMC, UNKNOWN, judge

CASES = {
    "Cycling NK cells": "NK Proliferating",
    "Proliferating NK cell": "NK Proliferating",
    "proliferating CD8+ T cells": "CD8 Proliferating",
    "NK cell": "NK",
    "CD56bright NK cell": "NK_CD56bright",
    "CD8+ effector memory T cell": "CD8 TEM",
    "T cell": "T cell (CD4/CD8/other unclear)",
    "Cytotoxic CD4+ T cell": "CD4 CTL",
    "gamma-delta T cell": "gdT",
    "Classical monocyte": "CD14 Mono",
    "Neutrophil": NOT_PBMC,
    "Regulatory T cell": "Treg",
    "naive CD4+ T cell": "CD4 Naive",
    "memory B cell": "B memory",
    "B cell": "B cell (subtype unclear)",
    "Monocyte": "Monocyte (subtype unclear)",
    "Reticulocyte": "Eryth",
    "AXL+ SIGLEC6+ DC": "ASDC",
    "Unknown": UNKNOWN,
    "Plasma cell": "Plasmablast",
    # found too lenient in the first spot-check:
    "Dendritic cells (myeloid/conventional DCs)": "Dendritic cell (subtype unclear)",
    "Proliferating cells (likely cycling lymphocytes/progenitor cells)": "Proliferating lymphocyte (lineage unclear)",
    "Proliferating lymphocytes (cycling T/NK cells)": "Proliferating lymphocyte (lineage unclear)",
    "Proliferating T cells (cycling lymphocytes)": "Proliferating lymphocyte (lineage unclear)",
    "Proliferating NK cells (cycling)": "NK Proliferating",
    "CD8+ cytotoxic T cells (or NK-like/NKT cells)": "CD8 T cell (subtype unclear)",
    "CD8+ cytotoxic T cell / NK-like T cell": "CD8 T cell (subtype unclear)",
    "cDC2 (CD1C+ dendritic cell)": "cDC2",
    "Naive T cell (likely naive CD4+ T cell)": "CD4 Naive",
    "Lymphocyte (likely naive T cell)": "T cell (CD4/CD8/other unclear)",
    # the 3 disagreements from the 30-row hand-check:
    "CD8+ memory/effector T cells (KLRB1+ GZMK+ subset, possibly MAIT-like)": "CD8 TEM",
    "cytotoxic CD8+ T cell": "CD8 T cell (subtype unclear)",
    "CD8+ cytotoxic T cell": "CD8 T cell (subtype unclear)",
    # rule 2: effector wording still maps to TEM
    "CD8+ effector T cells": "CD8 TEM",
    "CD8+ effector memory T cell": "CD8 TEM",
}

misses = 0
for answer, expected in CASES.items():
    got = judge(answer)
    if got != expected:
        misses += 1
        print(f"MISS  {answer!r}: expected {expected!r}, got {got!r}")
print(f"{len(CASES) - misses}/{len(CASES)} correct")
