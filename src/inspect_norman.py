"""Summarize Norman et al. 2019 (scPerturb NormanWeissman2019_filtered.h5ad) without loading the matrix.

Opens the file in backed (read-only, on-disk) mode, so only cell metadata goes into memory.
Prints: cells x genes, metadata columns, control cells, and perturbations split into single-gene vs
two-gene (pair) perturbations, with cells per perturbation.

Usage: python src/inspect_norman.py
"""
import anndata as ad
import pandas as pd

PATH = "data/norman2019/NormanWeissman2019_filtered.h5ad"

a = ad.read_h5ad(PATH, backed="r")
obs = a.obs
print(f"{a.n_obs:,} cells x {a.n_vars:,} genes | X stored as {type(a.X).__name__} (backed, not loaded)")
print("\nobs columns:", list(obs.columns))
print("var columns:", list(a.var.columns))
for col in ["perturbation", "nperts", "perturbation_type", "disease", "cell_line", "celltype", "organism"]:
    if col in obs:
        vc = obs[col].value_counts(dropna=False)
        print(f"\n{col}: {len(vc)} distinct" + ("" if len(vc) > 12 else f" -> {vc.to_dict()}"))

# ---- perturbation summary --------------------------------------------------------------------
# Labels: "control" (nperts 0), single gene "KLF1" (nperts 1), pair "SET_KLF1" (nperts 2).
# Note: Norman et al. 2019 is a CRISPR *activation* (CRISPRa) screen - perturbed genes are switched on.
per = obs.groupby("perturbation", observed=True).agg(cells=("nperts", "size"), nperts=("nperts", "first"))
kind = per.nperts.map({0: "control", 1: "single", 2: "pair"})
print("\n=== SUMMARY ===")
print(f"cells: {a.n_obs:,} | genes measured: {a.n_vars:,} | cell line: K562 (CRISPRa screen)")
for k in ["control", "single", "pair"]:
    sub = per[kind == k]
    print(f"{k:>8}: {len(sub):>3} perturbation label(s), {sub.cells.sum():>6,} cells"
          + (f" | cells per perturbation: median {sub.cells.median():.0f}, min {sub.cells.min()}, max {sub.cells.max()}"
             if k != "control" else ""))
singles = set(per.index[kind == "single"])
pairs = [p.split("_") for p in per.index[kind == "pair"]]
both = sum(all(g in singles for g in p) for p in pairs)
print(f"pairs whose two genes are BOTH also perturbed alone: {both}/{len(pairs)}")
print(f"distinct genes perturbed (singles or in pairs): {len(singles | {g for p in pairs for g in p})}")
print(f"perturbations with < 50 cells: {(per[kind != 'control'].cells < 50).sum()}")
