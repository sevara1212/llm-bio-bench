"""Step 5: map free-text predictions to canonical cell types and score them.

Usage: python src/score.py [pbmc|hao]   (default pbmc) -> results/scores_<dataset>.csv

pbmc: regex RULES below map answers to the 8 PBMC3k types.
hao:  the LLM judge (judge.py, run it first) maps answers to the 30 Hao types.

Each answer gets a score:
  1   = correct cell type
  0.5 = right lineage only - either too vague ("T cell" for CD4 T cells) or,
        on hao, the wrong subtype of the right lineage (CD4 TCM for CD4 Naive)
  0   = wrong
strict = share scoring 1; lenient = mean score.
"""
import glob
import json
import re
import sys

import pandas as pd

DATASET = sys.argv[1] if len(sys.argv) > 1 else "pbmc"

# --- pbmc: regex mapping -----------------------------------------------------
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


def normalize_pbmc(pred):
    if not isinstance(pred, str):
        return "no answer"
    p = pred.lower()
    for label, pattern in RULES:
        if re.search(pattern, p):
            return label
    return "other"


def score_pbmc(row):
    if row.normalized == row.answer:
        return 1.0
    return 0.5 if row.answer in PARTIAL.get(row.normalized, ()) else 0.0


# --- hao: judge mapping + lineage (celltype.l1) partial credit ---------------
def hao_scorer():
    from judge import CACHE, VAGUE

    cache = json.load(open(CACHE))
    l1 = pd.read_csv("data/hao_labels.csv").set_index("cell_type").l1.to_dict()
    no_partial = {"other"}  # Eryth / HSPC / ILC / Platelet share an l1 but aren't one lineage

    def lineage(label):
        if label in l1:
            return l1[label]
        return VAGUE.get(label)  # e.g. "CD4 T cell (subtype unclear)" -> "CD4 T"; None for unknown

    def normalize(pred):
        if not isinstance(pred, str):
            return "no answer"
        return cache.get(pred, "NOT JUDGED")

    def score(row):
        if row.normalized == row.answer:
            return 1.0
        truth, pred = l1[row.answer], lineage(row.normalized)
        if pred is None or (truth in no_partial and pred != "Proliferating"):
            return 0.0
        if pred == truth or (pred == "T" and truth in {"CD4 T", "CD8 T", "other T"}):
            return 0.5
        if pred == "Proliferating" and row.answer.endswith("Proliferating"):
            return 0.5  # spotted the cycling state, hedged on lineage
        return 0.0

    return normalize, score, set(VAGUE)


if DATASET == "pbmc":
    normalize, score, vague_labels = normalize_pbmc, score_pbmc, set(PARTIAL)
else:
    normalize, score, vague_labels = hao_scorer()

frames = []
for path in sorted(glob.glob(f"results/*_{DATASET}.csv")):
    if path.split("/")[-1].startswith("scores_"):
        continue
    df = pd.read_csv(path)
    if "knob" not in df.columns:  # skip old-format files
        continue
    if "run" not in df.columns:
        df["run"] = 0
    df.insert(0, "model", path.split("/")[-1].removesuffix(f"_{DATASET}.csv"))
    frames.append(df)
scores = pd.concat(frames, ignore_index=True)
scores["normalized"] = scores.predicted.map(normalize)
if (scores.normalized == "NOT JUDGED").any():
    sys.exit(f"{(scores.normalized == 'NOT JUDGED').sum()} answers not judged yet: run python src/judge.py")
scores["score"] = scores.apply(score, axis=1)
scores["strict"] = (scores.score == 1).astype(float)
scores["vague"] = scores.normalized.isin(vague_labels).astype(float)
scores.to_csv(f"results/scores_{DATASET}.csv", index=False)

pd.set_option("display.width", 160)
print(f"[{DATASET}] runs per question: {scores.groupby('model').run.nunique().to_dict()}\n")
print("Overall accuracy (strict = exact cell type; lenient = 0.5 credit for right lineage)\n")
print(scores.groupby(["model", "format"])[["strict", "score"]].mean().round(3)
      .rename(columns={"score": "lenient"}).unstack().to_string())
print("\nVague-answer rate (named a lineage but no subtype)\n")
print(scores.groupby("model").vague.mean().round(3).to_string())
print("\nCost and tokens per answer\n")
print(scores.groupby("model").agg(answers=("id", "size"), total_cost=("cost_usd", "sum"),
                                  cost_per_answer=("cost_usd", "mean"),
                                  output_tokens=("completion_tokens", "mean")).round(4).to_string())
if DATASET == "pbmc":
    print("\nUnmatched predictions (check whether RULES should cover them):")
    print(scores[scores.normalized == "other"].predicted.value_counts().head(15).to_string())
