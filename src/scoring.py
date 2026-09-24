"""Answer -> label mapping and per-answer scores, shared by score.py and compare_baselines.py.

A score function takes a row with .answer (true label) and .normalized (mapped prediction):
1 = correct, 0.5 = right lineage only, 0 = wrong.
"""
import json
import re

import pandas as pd

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
    from judge import CACHE, VAGUE, rule_match

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
        return rule_match(pred) or cache.get(pred, "NOT JUDGED")

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


def scorer(dataset):
    """Return (normalize, score, vague_labels) for 'pbmc' or 'hao'."""
    if dataset == "pbmc":
        return normalize_pbmc, score_pbmc, set(PARTIAL)
    return hao_scorer()
