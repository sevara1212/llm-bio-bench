"""Task 2 (Norman et al. 2019, CRISPRa in K562): single-gene perturbation -> target gene direction.

Usage: python src/task2_build.py

1. Pseudobulk: mean log-normalised expression (counts / cell total x 1e4, then natural log1p) per
   single-gene perturbation and for control. The matrix is read in blocks of genes (it is stored by
   gene), so only one block of genes is in memory at a time. Pair-perturbation cells are not used.
2. Differential expression for every (activated gene X, target gene Y):
     log2FC = (mean log-norm in X cells - mean in control cells) / ln 2   (difference of means, log2 units)
     Welch's t-test on per-cell log-normalised values, computed exactly from per-group sums and sums
     of squares; Benjamini-Hochberg across the tested Y genes of each X (one family per perturbation).
   Tested Y: expressed (count > 0) in >= 10% of control cells; Y != X.
3. Labels: up = log2FC > 0.5 & padj < 0.05; down = log2FC < -0.5 & padj < 0.05;
           no_change = |log2FC| < 0.1 & padj > 0.5; everything else is dropped.
4. Sample ~200 questions (seed 0): balanced over up / down / no_change, round-robin over X genes,
   at most 5 questions per X and at most 3 per target gene Y; Y with clone-style names (RP11-...,
   AC012345.1, ...) are excluded. up and down are sampled first; each no_change question is then
   matched on expression to one up/down question: its Y must be detected in a similar share of
   control cells (within 0.02, widened to 0.05 if needed), so "rarely detected -> no_change" is not a
   shortcut (without matching, no_change Y had a median detection of 32% vs 94-96% for up/down).

Outputs:
  data/norman2019/task2_pseudobulk_singles.csv.gz   mean log-norm expression, groups x all genes (gitignored)
  data/norman2019/task2_de_singles.csv.gz            every tested (X, Y) with log2FC, t, p, padj, label (gitignored)
  data/norman2019/task2_ctrl_corr.npz                control-cell Pearson r of each X with every gene (gitignored;
                                                     used by the co-expression baseline and agent tool)
  data/task2_questions.json                          the sampled questions
"""
import json
import math
import random

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

PATH = "data/norman2019/NormanWeissman2019_filtered.h5ad"
BLOCK = 2000
MIN_CTRL_FRAC = 0.10
N_QUESTIONS, MAX_PER_X, MAX_PER_Y, SEED = 200, 5, 3, 0
CLONE_NAME = r"^(RP\d+-|CT[ABCD]-|CTD-|AC\d{6}|AL\d{6}|AP\d{6}|AF\d{6}|XXbac-|KB-|LA16c-)"
MATCH_TOL = (0.02, 0.05)
# Norman's labels use old names for three genes; questions use current HGNC symbols.
CURRENT_SYMBOL = {"C19orf26": "CBARP", "C3orf72": "FOXL2NB", "KIAA1804": "MAP3K21"}
QUESTION = ("In K562 cells (a human leukemia cell line), gene {X} is activated with CRISPR activation "
            "(CRISPRa). What happens to the expression of gene {Y}? Answer up, down, or no_change, "
            'as JSON with a confidence 0-1.')

a = ad.read_h5ad(PATH, backed="r")
obs = a.obs
cells = np.flatnonzero((obs.nperts <= 1).to_numpy())  # control + single-gene perturbations
labels = obs.perturbation.astype(str).to_numpy()[cells]
groups = ["control"] + sorted(set(labels) - {"control"})
gidx = {g: i for i, g in enumerate(groups)}
cell_group = np.array([gidx[l] for l in labels])
n_cells, n_genes = len(cells), a.n_vars
G = sp.csr_matrix((np.ones(n_cells), (cell_group, np.arange(n_cells))), shape=(len(groups), n_cells))
n_per_group = np.asarray(G.sum(axis=1)).ravel()
ctrl = cell_group == 0
genes = a.var_names.to_numpy()
x_names = [CURRENT_SYMBOL.get(g, g) for g in groups[1:]]  # measured-gene name of each activated gene
x_cols = [int(np.flatnonzero(genes == x)[0]) if x in set(genes) else None for x in x_names]
print(f"{n_cells:,} cells ({ctrl.sum():,} control), {len(groups) - 1} single perturbations, {n_genes:,} genes")


def blocks():
    for start in range(0, n_genes, BLOCK):
        stop = min(start + BLOCK, n_genes)
        yield start, stop, sp.csc_matrix(a.X[:, start:stop])[cells]


# ---- pass 1: cell totals (check against obs.ncounts) and raw counts of the activated genes --------
totals = np.zeros(n_cells)
x_raw = np.zeros((n_cells, len(x_cols)))
for start, stop, M in blocks():
    totals += np.asarray(M.sum(axis=1)).ravel()
    for j, c in enumerate(x_cols):
        if c is not None and start <= c < stop:
            x_raw[:, j] = M[:, c - start].toarray().ravel()
stored = obs.ncounts.to_numpy()[cells]
print(f"cell totals match obs.ncounts: {np.allclose(totals, stored)} (max abs diff {np.abs(totals - stored).max():.1f})")
scale = 1e4 / totals
x_log = np.log1p(x_raw * scale[:, None])  # activated genes, log-normalised, all used cells
xc = x_log[ctrl]  # control cells only, for co-expression

# ---- pass 2: per-group sums, sums of squares, control detection, control cross-products ----------
S1 = np.zeros((len(groups), n_genes))
S2 = np.zeros((len(groups), n_genes))
ctrl_frac = np.zeros(n_genes)
cross = np.zeros((len(x_cols), n_genes))  # sum over control cells of x_X * y_Y
for start, stop, M in blocks():
    L = sp.diags(scale) @ M.tocsr()
    L.data = np.log1p(L.data)
    S1[:, start:stop] = (G @ L).toarray()
    S2[:, start:stop] = (G @ L.multiply(L)).toarray()
    Lc = L[ctrl]
    ctrl_frac[start:stop] = np.asarray((Lc > 0).sum(axis=0)).ravel() / ctrl.sum()
    cross[:, start:stop] = np.asarray(Lc.T @ xc).T
mean = S1 / n_per_group[:, None]
var = np.maximum(S2 / n_per_group[:, None] - mean ** 2, 0) * n_per_group[:, None] / (n_per_group[:, None] - 1)
pd.DataFrame(mean, index=groups, columns=genes).to_csv("data/norman2019/task2_pseudobulk_singles.csv.gz")

# control-cell Pearson r between each activated gene X and every gene Y
nc = ctrl.sum()
mx, sx = xc.mean(axis=0), xc.std(axis=0)
my, sy = mean[0], np.sqrt(var[0] * (nc - 1) / nc)
with np.errstate(divide="ignore", invalid="ignore"):
    corr = (cross / nc - np.outer(mx, my)) / np.outer(sx, sy)
np.savez_compressed("data/norman2019/task2_ctrl_corr.npz", corr=corr.astype(np.float32), x_genes=np.array(x_names),
                    genes=genes, ctrl_frac=ctrl_frac.astype(np.float32),
                    x_ctrl_frac=(xc > 0).mean(axis=0).astype(np.float32))

# ---- differential expression, labels --------------------------------------------------------------
tested = ctrl_frac >= MIN_CTRL_FRAC
rows = []
for j, (label, x) in enumerate(zip(groups[1:], x_names)):
    g = j + 1
    keep = tested & (genes != x)
    d = mean[g, keep] - mean[0, keep]
    va, vb, na, nb = var[g, keep] / n_per_group[g], var[0, keep] / nc, n_per_group[g], nc
    se = np.sqrt(va + vb)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, d / se, 0.0)
        df = np.where(se > 0, (va + vb) ** 2 / (va ** 2 / (na - 1) + vb ** 2 / (nb - 1)), 1.0)
    p = np.where(se > 0, 2 * stats.t.sf(np.abs(t), df), 1.0)
    padj = stats.false_discovery_control(p, method="bh")
    rows.append(pd.DataFrame({
        "X": x, "X_norman_label": label, "Y": genes[keep], "log2FC": d / math.log(2), "t": t, "p": p,
        "padj": padj, "ctrl_frac_Y": ctrl_frac[keep], "n_X_cells": int(na), "n_ctrl_cells": int(nb),
        "coexpr_r_ctrl": corr[j, keep]}))
de = pd.concat(rows, ignore_index=True)
de["label"] = np.select(
    [(de.log2FC > 0.5) & (de.padj < 0.05), (de.log2FC < -0.5) & (de.padj < 0.05),
     (de.log2FC.abs() < 0.1) & (de.padj > 0.5)], ["up", "down", "no_change"], default="dropped")
de.to_csv("data/norman2019/task2_de_singles.csv.gz", index=False)

# CRISPRa sanity check: the activated gene itself should go up in its own perturbation.
self_fc = pd.Series({x: (mean[j + 1, c] - mean[0, c]) / math.log(2)
                     for j, (x, c) in enumerate(zip(x_names, x_cols)) if c is not None})
print(f"\nCRISPRa check - activated gene's own log2FC: median {self_fc.median():.2f}, "
      f"> 0.5 in {(self_fc > 0.5).sum()}/{len(self_fc)} perturbations (MAP3K21 not measured)")
print(f"tested target genes (>= {MIN_CTRL_FRAC:.0%} of control cells): {tested.sum():,}")
print(f"(X, Y) pairs tested: {len(de):,}\nlabels over all pairs:")
print(de.label.value_counts().to_string())

# ---- sample the questions ------------------------------------------------------------------------
rng = random.Random(SEED)
cand = de[(de.label != "dropped") & ~de.Y.str.match(CLONE_NAME)]
print(f"\nexcluded clone-style Y names: {(de.label != 'dropped').sum() - len(cand):,} labelled pairs")
per_x, per_y, chosen = {}, {}, []


def take(i):
    chosen.append(i)
    per_x[de.at[i, "X"]] = per_x.get(de.at[i, "X"], 0) + 1
    per_y[de.at[i, "Y"]] = per_y.get(de.at[i, "Y"], 0) + 1


def allowed(i):
    return per_x.get(de.at[i, "X"], 0) < MAX_PER_X and per_y.get(de.at[i, "Y"], 0) < MAX_PER_Y


# up and down: one question per class per round, rotating over X genes (X moves to the back when used)
quota = {"up": N_QUESTIONS // 3 + (N_QUESTIONS % 3 > 0), "down": N_QUESTIONS // 3 + (N_QUESTIONS % 3 > 1)}
pool = {lab: {x: rng.sample(list(sub.index), len(sub)) for x, sub in cand[cand.label == lab].groupby("X")}
        for lab in quota}
xs = sorted(cand.X.unique())
rng.shuffle(xs)
have = {lab: 0 for lab in quota}
progress = True
while progress and any(have[l] < quota[l] for l in quota):
    progress = False
    for lab in quota:
        if have[lab] >= quota[lab]:
            continue
        for x in xs:
            idx = next((i for i in pool[lab].get(x, []) if allowed(i)), None)
            if idx is not None and per_x.get(x, 0) < MAX_PER_X:
                pool[lab][x].remove(idx)
                take(idx)
                have[lab] += 1
                progress = True
                xs.append(xs.pop(xs.index(x)))
                break

# no_change: match each one to a random up/down question on Y's control detection rate
n_nc = N_QUESTIONS - len(chosen)
anchors = rng.sample(chosen, n_nc)
nc_pool = cand[cand.label == "no_change"]
unmatched = 0
for a_i in anchors:
    target = de.at[a_i, "ctrl_frac_Y"]
    pick = None
    for tol in MATCH_TOL:
        near = nc_pool[(nc_pool.ctrl_frac_Y - target).abs() <= tol]
        near = near[[allowed(i) for i in near.index]]
        if len(near):
            fewest = near.X.map(lambda x: per_x.get(x, 0)).min()  # spread over X genes
            near = near[near.X.map(lambda x: per_x.get(x, 0)) == fewest]
            pick = rng.choice(list(near.index))
            break
    if pick is None:
        unmatched += 1
        continue
    take(pick)
    nc_pool = nc_pool.drop(pick)
print(f"no_change questions without an expression match within {MATCH_TOL[-1]}: {unmatched}")

q = de.loc[chosen].reset_index(drop=True)
questions = [{"id": f"t2_{i:03d}", "X": r.X, "X_norman_label": r.X_norman_label, "Y": r.Y, "label": r.label,
              "log2FC": round(float(r.log2FC), 4), "padj": float(r.padj), "ctrl_frac_Y": round(float(r.ctrl_frac_Y), 4),
              "coexpr_r_ctrl": None if np.isnan(r.coexpr_r_ctrl) else round(float(r.coexpr_r_ctrl), 5),
              "question": QUESTION.format(X=r.X, Y=r.Y)} for i, r in enumerate(q.itertuples())]
json.dump(questions, open("data/task2_questions.json", "w"), indent=1)
print(f"\nSaved {len(questions)} questions to data/task2_questions.json")
print("label counts:", q.label.value_counts().to_dict())
print(f"distinct X genes: {q.X.nunique()} | questions per X: max {q.X.value_counts().max()}, "
      f"median {q.X.value_counts().median():.0f}")
print(f"distinct Y genes: {q.Y.nunique()} | questions per Y: max {q.Y.value_counts().max()}")
