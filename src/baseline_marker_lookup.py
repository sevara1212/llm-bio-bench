"""Marker-lookup baseline (no LLM): pick the cell type whose known markers overlap most with the genes.

Usage: python src/baseline_marker_lookup.py <pbmc|hao> <panglaodb|cellmarker>
Writes results/<db>-lookup_<dataset>.csv in the same format as run_eval.py, so score.py scores it
exactly like a model.

Method, identical for both databases:
  1. Candidate cell types = database types mapped to this project's labels (tables below /
     configs/cellmarker_to_hao.csv). Unmapped types are not candidates.
  2. Known markers = each type's top TOP_N genes by evidence (CellMarker: number of distinct papers;
     PanglaoDB: human sensitivity), so big generic lists (B cell: 1,298 genes) don't win by size.
  3. Ensembl IDs are converted to symbols with the dataset's own ID table (a lookup tool doesn't need
     to have memorized IDs).
  4. Pick the type with the most given genes among its known markers. Ties -> the type whose matched
     markers' ranks sum lowest (stronger markers win), then alphabetical. No overlap -> Unknown.
     The labels of all types tied on the top count are saved too (tied_labels), so the baseline can
     also be scored with fractional credit: 1/k for each of k tied labels (compare_baselines.py).

CellMarker rows used: human, tissue class Blood, disease Normal, EXCLUDING Hao et al. 2021
(PMID 34062119), whose own marker lists are in CellMarker and would leak the answer key.
"""
import csv
import json
import os
import sys

import pandas as pd

from judge import NOT_PBMC, UNKNOWN

DATASET, DB = (sys.argv[1:3] + ["hao", "cellmarker"][len(sys.argv[1:3]):])[:2]
TOP_N = 50
HAO_PMID = 34062119

PANGLAO = "data/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
CELLMARKER_RAW = "data/markers/human_cell_marker.txt"
CELLMARKER_BLOOD = "data/markers/cellmarker_blood_normal.csv"  # filtered extract, made on first run

U = {"B": "B cell (subtype unclear)", "CD4": "CD4 T cell (subtype unclear)", "CD8": "CD8 T cell (subtype unclear)",
     "T": "T cell (CD4/CD8/other unclear)", "Mono": "Monocyte (subtype unclear)", "DC": "Dendritic cell (subtype unclear)"}
# PanglaoDB type -> Hao label. Types PanglaoDB can't split map to the lineage-only label (lenient credit).
PANGLAO_TO_HAO = {
    "Platelets": "Platelet", "Megakaryocytes": "Platelet", "Plasmacytoid dendritic cells": "pDC",
    "Dendritic cells": U["DC"], "Monocytes": U["Mono"], "B cells": U["B"], "B cells naive": "B naive",
    "B cells memory": "B memory", "Plasma cells": "Plasmablast", "T cells": U["T"], "T memory cells": U["T"],
    "Natural killer T cells": U["T"], "T helper cells": U["CD4"], "T follicular helper cells": U["CD4"],
    "T regulatory cells": "Treg", "T cytotoxic cells": U["CD8"], "Gamma delta T cells": "gdT", "NK cells": "NK",
    "Nuocytes": "ILC", "Erythroid-like and erythroid precursor cells": "Eryth", "Erythroblasts": "Eryth",
    "Reticulocytes": "Eryth", "Hematopoietic stem cells": "HSPC",
    **{t: NOT_PBMC for t in ["Macrophages", "Red pulp macrophages", "Mast cells", "Basophils", "Neutrophils",
                             "Eosinophils", "Myeloid-derived suppressor cells"]},
}
# Hao label -> PBMC3k label, worded so score.py's PBMC regex RULES read it as intended.
HAO_TO_PBMC = {
    **{l: "CD4 T cells" for l in ["CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating", "Treg", U["CD4"]]},
    **{l: "CD8 T cells" for l in ["CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", U["CD8"]]},
    **{l: "T cells (unspecified)" for l in ["gdT", "MAIT", "dnT", U["T"]]},
    **{l: "NK cells" for l in ["NK", "NK_CD56bright", "NK Proliferating"]},
    **{l: "B cells" for l in ["B naive", "B memory", "B intermediate", "Plasmablast", U["B"]]},
    **{l: "Dendritic cells" for l in ["cDC1", "cDC2", "pDC", "ASDC", U["DC"]]},
    "CD14 Mono": "CD14 Monocytes", "CD16 Mono": "FCGR3A Monocytes", U["Mono"]: "Monocytes (unspecified)",
    "Platelet": "Megakaryocytes",
}


def load_markers():
    """Return ({db cell type: [top-N marker symbols, strongest first]}, {db cell type: Hao label})."""
    if DB == "panglaodb":
        db = pd.read_csv(PANGLAO, sep="\t")
        db = db[db.species.str.contains("Hs") & db["cell type"].isin(PANGLAO_TO_HAO)]
        db = db.assign(symbol=db["official gene symbol"].str.upper())
        db = db.sort_values(["cell type", "sensitivity_human", "canonical marker", "symbol"],
                            ascending=[True, False, False, True])
        markers = {ct: g.drop_duplicates().head(TOP_N).tolist() for ct, g in db.groupby("cell type").symbol}
        return markers, PANGLAO_TO_HAO

    if not os.path.exists(CELLMARKER_BLOOD):
        cols = ["species", "tissue_class", "disease", "cell_name", "symbol", "pmid"]
        parts = [ch[(ch.species == "Human") & (ch.tissue_class == "Blood") & (ch.disease == "Normal")]
                 for ch in pd.read_csv(CELLMARKER_RAW, sep="\t", usecols=cols, chunksize=500_000, low_memory=False)]
        pd.concat(parts).to_csv(CELLMARKER_BLOOD, index=False)
    db = pd.read_csv(CELLMARKER_BLOOD)
    db = db[(db.pmid != HAO_PMID) & db.symbol.notna()]
    mapping = dict(csv.reader(open("configs/cellmarker_to_hao.csv")))
    mapping.pop("cellmarker_cell_name", None)
    db = db[db.cell_name.isin(mapping)].assign(symbol=lambda d: d.symbol.str.upper())
    support = (db.groupby(["cell_name", "symbol"]).agg(papers=("pmid", "nunique"), rows=("pmid", "size"))
               .reset_index().sort_values(["cell_name", "papers", "rows", "symbol"], ascending=[True, False, False, True]))
    markers = {ct: g.head(TOP_N).tolist() for ct, g in support.groupby("cell_name").symbol}
    return markers, mapping



def lookup(genes, markers, ensembl_to_symbol):
    """Return (picked database type or None, overlap count, all database types tied at that count).

    Pick = most given genes among the type's markers; ties broken by the lowest sum of the matched
    markers' ranks (stronger markers win), then alphabetically.
    """
    symbols = {str(ensembl_to_symbol.get(g, g)).upper() for g in genes}
    hits = {ct: [i for i, m in enumerate(ranked) if m in symbols] for ct, ranked in markers.items()}
    top = max(len(h) for h in hits.values())
    if top == 0:
        return None, 0, []
    tied = sorted(ct for ct, h in hits.items() if len(h) == top)
    pick = min(tied, key=lambda ct: (sum(hits[ct]), ct))
    return pick, top, tied


def main():
    markers, to_hao = load_markers()
    id_tables = [pd.read_csv(f"data/{DATASET}_markers.csv"), pd.read_csv(f"data/{DATASET}_background_genes.csv")]
    ensembl_to_symbol = {e: g for t in id_tables for e, g in zip(t.ensembl, t.gene)}

    def to_label(ct):
        if ct is None:
            return UNKNOWN if DATASET == "hao" else "no marker overlap"
        hao = to_hao[ct]
        return hao if DATASET == "hao" else HAO_TO_PBMC.get(hao, "other")

    questions = json.load(open(f"data/questions_{DATASET}.json"))
    # Match the models' runs so score.py's same-questions filter keeps every run (the baseline is
    # deterministic, so all runs are identical).
    runs = sorted(pd.read_csv(f"results/claude-sonnet-5_{DATASET}.csv").run.unique())
    rows = []
    for q in questions:
        ct, overlap, tied = lookup(q["genes"], markers, ensembl_to_symbol)
        # Distinct labels among tied types, for fractional credit (1/k each). Several database
        # names can map to one label (e.g. 7 CellMarker NK names), so count labels, not names.
        tied_labels = sorted({to_label(t) for t in tied}) if tied else [to_label(None)]
        for run in runs:
            rows.append({
                "id": q["id"], "format": q["format"], "knob": q["knob"], "level": q["level"],
                "replicate": q["replicate"], "run": run, "provider": "none", "answer": q["answer"],
                "predicted": to_label(ct), "confidence": None, "raw": ct or "no overlap",
                "finish_reason": None, "prompt_tokens": 0, "completion_tokens": 0, "reasoning_tokens": 0,
                "cost_usd": 0.0, "marker_overlap": overlap, "n_tied_labels": len(tied_labels),
                "tied_labels": json.dumps(tied_labels),
            })
    out = pd.DataFrame(rows)
    path = f"results/{DB}-lookup_{DATASET}.csv"
    out.to_csv(path, index=False)
    one = out[out.run == runs[0]]
    sizes = pd.Series({ct: len(m) for ct, m in markers.items()})
    print(f"[{DB} / {DATASET}] {len(markers)} candidate cell types, markers per type: median {sizes.median():.0f}, "
          f"min {sizes.min()} ({sizes.idxmin()})")
    print(f"Wrote {len(out)} rows ({len(one)} questions x {len(runs)} run(s)) to {path}")
    print(f"no marker overlap: {(one.marker_overlap == 0).mean():.1%} | top count tied between 2+ labels: "
          f"{(one.n_tied_labels > 1).mean():.1%} (median k when tied: "
          f"{one[one.n_tied_labels > 1].n_tied_labels.median():.0f})")


if __name__ == "__main__":
    main()
