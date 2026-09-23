"""Step 5: map free-text predictions to canonical cell types and score them.

Each answer gets a score: 1 = correct, 0.5 = right lineage but too vague
(e.g. "T cell" for CD4 T cells, "monocyte" for CD14 Monocytes), 0 = wrong.
Writes results/scores.csv (one row per answer) and prints accuracy tables.
"""
import glob
import re

import pandas as pd

# Checked in order; first match wins. Order matters: "non-classical" before
# "classical", CD8 before NK (CD8 answers often say "NK-like"), specific T before generic T.
RULES = [
    ("Megakaryocytes", r"megakaryocyte|platelet"),
    ("B cells", r"\bb[ -]?(cell|lymph)|plasma ?cell"),
    ("FCGR3A Monocytes", r"fcgr3a|cd16|non-?classical|nonclassical"),
    ("CD14 Monocytes", r"cd14|classical mono"),
    ("Dendritic cells", r"dendritic|\bc?dc\d?\b"),
    ("CD8 T cells", r"cd8|cytotoxic t"),
    ("CD4 T cells", r"cd4|helper t|\bth\d"),
    ("NK cells", r"\bnk\b|natural killer"),
    ("T cells (unspecified)", r"\bt[ -]?(cell|lymph)"),
    ("Monocytes (unspecified)", r"monocyte|myeloid"),
]
PARTIAL = {
    "T cells (unspecified)": {"CD4 T cells", "CD8 T cells"},
    "Monocytes (unspecified)": {"CD14 Monocytes", "FCGR3A Monocytes"},
}


def normalize(pred):
    if not isinstance(pred, str):
        return "no answer"
    p = pred.lower()
    for label, pattern in RULES:
        if re.search(pattern, p):
            return label
    return "other"


def score(row):
    if row.normalized == row.answer:
        return 1.0
    return 0.5 if row.answer in PARTIAL.get(row.normalized, ()) else 0.0


frames = []
for path in sorted(glob.glob("results/*_pbmc.csv")):
    df = pd.read_csv(path)
    if "knob" not in df.columns:  # skip old-format files
        continue
    df.insert(0, "model", path.split("/")[-1].removesuffix("_pbmc.csv"))
    frames.append(df)
scores = pd.concat(frames, ignore_index=True)
scores["normalized"] = scores.predicted.map(normalize)
scores["score"] = scores.apply(score, axis=1)
scores["strict"] = (scores.score == 1).astype(float)
scores.to_csv("results/scores.csv", index=False)

pd.set_option("display.width", 160)
print("Overall accuracy (strict = exact cell type; lenient = partial credit)\n")
print(scores.groupby(["model", "format"])[["strict", "score"]].mean().round(3).unstack().to_string())
print("\nCost and tokens per answer\n")
print(scores.groupby("model").agg(answers=("id", "size"), total_cost=("cost_usd", "sum"),
                                  cost_per_answer=("cost_usd", "mean"),
                                  output_tokens=("completion_tokens", "mean")).round(4).to_string())
print("\nUnmatched predictions (check whether RULES should cover them):")
print(scores[scores.normalized.isin(["other"])].predicted.value_counts().head(15).to_string())
