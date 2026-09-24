"""Baselines vs LLMs (run score.py hao --skip-knob noise and score.py pbmc first).

Prints: overall table incl. the random-guess baseline; lookup baselines scored two ways (tie-break
pick vs fractional credit 1/k over the k labels tied on the top count); the comparison restricted to
the cell types each database can name; strict accuracy per difficulty knob; per lineage.
"""
import json
from types import SimpleNamespace

import pandas as pd

from scoring import scorer

s = pd.read_csv("results/scores_hao.csv", dtype={"level": str})
labels = pd.read_csv("data/hao_labels.csv").set_index("cell_type")
L1 = labels.l1
LOOKUPS = {"panglaodb-lookup": "PanglaoDB", "cellmarker-lookup": "CellMarker"}
pd.set_option("display.width", 200)

# Random guess over the 30 fine types: strict 1/30; lenient adds 0.5 for a same-lineage wrong pick.
rand_lenient = sum(1 / 30 + 0.5 * sum(L1[o] == L1[t] and o != t and L1[t] != "other" for o in L1.index) / 30
                   for t in L1.index) / 30

overall = s.groupby(["model", "format"]).agg(strict=("strict", "mean"), lenient=("score", "mean")).unstack()
overall.columns = [f"{m} ({f})" for m, f in overall.columns]
overall.loc["random guess"] = [1 / 30, 1 / 30, rand_lenient, rand_lenient]
print(f"n = {s.groupby('model').size().iloc[0]} answers per model (300 per format)\n")
print("OVERALL (share of answers)")
print(overall.round(3).to_string())



def fractional(scores, dataset):
    """Lookup rows only: strict and lenient with 1/k credit spread over the k tied labels."""
    normalize, score, _ = scorer(dataset)
    rows = scores[scores.model.isin(LOOKUPS)].copy()
    strict, lenient = [], []
    for r in rows.itertuples():
        labs = [normalize(l) for l in json.loads(r.tied_labels)]
        strict.append(sum(l == r.answer for l in labs) / len(labs))
        lenient.append(sum(score(SimpleNamespace(answer=r.answer, normalized=l)) for l in labs) / len(labs))
    rows["strict_frac"], rows["lenient_frac"] = strict, lenient
    t = rows.groupby("model").agg(n=("id", "size"), strict_tiebreak=("strict", "mean"), strict_fractional=("strict_frac", "mean"),
                                  lenient_tiebreak=("score", "mean"), lenient_fractional=("lenient_frac", "mean"),
                                  share_tied=("n_tied_labels", lambda k: (k > 1).mean()))
    return t.round(3)


print("\nLOOKUP BASELINES: TIE-BREAK PICK vs FRACTIONAL CREDIT (1/k per tied label), Hao")
print(fractional(s, "hao").to_string())
p = pd.read_csv("results/scores_pbmc.csv", dtype={"level": str})
print("\nsame, PBMC3k (8 types):")
print(fractional(p, "pbmc").to_string())

# Fine types each database CAN name exactly, from its mapping table - fixed before seeing results.
# (Selecting types the lookup happened to predict would be post-hoc and favour the lookup.)
import csv  # noqa: E402

from baseline_marker_lookup import PANGLAO_TO_HAO  # noqa: E402
cellmarker_map = dict(csv.reader(open("configs/cellmarker_to_hao.csv")))
CAN_NAME = {"panglaodb-lookup": set(PANGLAO_TO_HAO.values()), "cellmarker-lookup": set(cellmarker_map.values())}
for model, name in LOOKUPS.items():
    covered = sorted(CAN_NAME[model] & set(labels.index))
    sub = s[s.answer.isin(covered)]
    t = sub.groupby(["model", "format"]).strict.mean().unstack().round(3)
    print(f"\nRESTRICTED TO THE {len(covered)} TYPES {name} CAN NAME (from its mapping; n = {len(sub) // s.model.nunique()} per model)")
    print(t.to_string())

print("\nSTRICT BY DIFFICULTY KNOB (ensembl | symbol)")
print(s.pivot_table(index=["knob", "level"], columns=["format", "model"], values="strict").round(2).to_string())

print("\nSTRICT BY LINEAGE (both formats)")
print(s.assign(l1=s.answer.map(L1)).pivot_table(index="l1", columns="model", values="strict").round(2).to_string())
