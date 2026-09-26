"""LLM judge for the Hao dataset: map each free-text answer to one of the 30 labels.

Two steps, applied identically to every model's answers:
  1. rule_match(): exact match (after lower-casing, trimming, dropping a plural "s") against the
     label names and a short list of unambiguous synonyms. No LLM involved.
  2. Only answers rule_match() can't place go to the LLM judge.

The judge sees only the answer text - never which model wrote it, the true label or the genes -
so it maps wording, it doesn't grade. Each distinct answer string is judged once and cached in
results/judge_cache_hao.json.

Usage: python src/judge.py           judge answers not yet in the cache
       python src/judge.py --fresh   ignore the cache and judge everything again in one pass
"""
import glob
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import anthropic
import pandas as pd
from dotenv import load_dotenv

load_dotenv()

JUDGE_MODEL = "claude-sonnet-5"  # knows immunology well enough to map "CD56bright NK" -> NK_CD56bright
CACHE = "results/judge_cache_hao.json"

FINE = sorted(pd.read_csv("data/hao_labels.csv").cell_type)
# Shown next to labels whose names alone are ambiguous.
DESCRIBE = {
    "NK": "ordinary / CD56dim NK cells - use for a plain 'NK cell' answer",
    "NK_CD56bright": "CD56bright NK cells only",
    "NK Proliferating": "cycling / proliferating / MKI67+ NK cells",
    "CD4 Proliferating": "cycling / proliferating CD4 T cells",
    "CD8 Proliferating": "cycling / proliferating CD8 T cells",
    "ASDC": "AXL+ SIGLEC6+ dendritic cells",
    "dnT": "double-negative T cells",
    "gdT": "gamma-delta T cells",
    "MAIT": "mucosal-associated invariant T cells",
    "HSPC": "haematopoietic stem / progenitor cells",
    "Eryth": "erythroid cells, reticulocytes",
    "Platelet": "platelets, megakaryocytes",
    "Plasmablast": "plasmablasts, plasma cells",
}
# Right lineage but no subtype: scored as partial credit (0.5).
VAGUE = {
    "B cell (subtype unclear)": "B",
    "CD4 T cell (subtype unclear)": "CD4 T",
    "CD8 T cell (subtype unclear)": "CD8 T",
    "T cell (CD4/CD8/other unclear)": "T",
    "Monocyte (subtype unclear)": "Mono",
    "Dendritic cell (subtype unclear)": "DC",
    "Proliferating lymphocyte (lineage unclear)": "Proliferating",
}
NOT_PBMC = "Not found in PBMCs (e.g. neutrophil, eosinophil, endothelial, epithelial)"
UNKNOWN = "Unknown / low-quality / no cell type given"
OPTIONS = FINE + list(VAGUE) + [NOT_PBMC, UNKNOWN]

# Unambiguous synonyms for rule_match(), written already normalized (lower case, singular).
SYNONYMS = {
    "nk cell": "NK", "natural killer cell": "NK", "natural killer (nk) cell": "NK",
    "cd56bright nk cell": "NK_CD56bright",
    "cd14+ monocyte": "CD14 Mono", "cd14 monocyte": "CD14 Mono", "classical monocyte": "CD14 Mono",
    "cd16+ monocyte": "CD16 Mono", "cd16 monocyte": "CD16 Mono", "non-classical monocyte": "CD16 Mono",
    "platelet": "Platelet", "megakaryocyte": "Platelet",
    "plasmablast": "Plasmablast", "plasma cell": "Plasmablast",
    "plasmacytoid dendritic cell": "pDC", "plasmacytoid dendritic cell (pdc)": "pDC",
    "naive b cell": "B naive", "memory b cell": "B memory",
    "regulatory t cell": "Treg", "regulatory t cell (treg)": "Treg",
    "gamma-delta t cell": "gdT", "mucosal-associated invariant t cell": "MAIT",
    "b cell": "B cell (subtype unclear)", "t cell": "T cell (CD4/CD8/other unclear)",
    "monocyte": "Monocyte (subtype unclear)", "dendritic cell": "Dendritic cell (subtype unclear)",
    "unknown": UNKNOWN,
}


def _normalize(text):
    t = re.sub(r"\s+", " ", text.strip().lower()).rstrip(".")
    return t[:-1] if len(t) > 3 and t.endswith("s") and not t.endswith("ss") else t


RULE_TABLE = {**{_normalize(o): o for o in OPTIONS}, **SYNONYMS}


def rule_match(answer):
    """Label for an answer that is exactly a label name or synonym, else None."""
    return RULE_TABLE.get(_normalize(answer)) if isinstance(answer, str) else None


INSTRUCTIONS = f"""A model was shown marker genes from a cluster of human PBMCs and named the cell type.
Map its answer to exactly one label from this list:

Fine cell types (use one of these whenever the answer is specific enough):
{chr(10).join('- ' + x + (f'  ({DESCRIBE[x]})' if x in DESCRIBE else '') for x in FINE)}

Only lineage given (use when the answer does not pin down a fine type above):
{chr(10).join('- ' + x for x in VAGUE)}

Other:
- {NOT_PBMC}
- {UNKNOWN}

Rules: TEM = effector memory, TCM = central memory, CTL = cytotoxic. "Classical monocyte" = CD14 Mono,
"non-classical" = CD16 Mono. Cycling / proliferating answers take priority over any other NK/CD4/CD8 subtype. If the answer names
one lineage, use it ("proliferating NK" -> NK Proliferating). If it hedges between lineages or names
none ("cycling T/NK cells", "proliferating lymphocytes") -> Proliferating lymphocyte (lineage unclear).
"Proliferating T cells" without CD4/CD8 -> Proliferating lymphocyte (lineage unclear).
CD8 TEM only if the answer says effector memory, effector or TEM. "Cytotoxic CD8 T cell" alone does not
specify a subtype -> CD8 T cell (subtype unclear).
"Conventional", "myeloid" or "classical" DC without saying cDC1/cDC2 (or CLEC9A/XCR1 vs CD1C/FCER1A)
-> Dendritic cell (subtype unclear).
Be strict: only pick a fine type when the answer itself says it. Never pick one because it is the most
common subtype.
Judge the PRIMARY answer only: the cell type the answer commits to. Ignore hedges and alternatives such as
"possibly X", "maybe X", "or Y", "X-like", "e.g. X" - never let a hedge pick the label. Words describing the
cell's state (proliferating, cycling, naive, memory, effector) are part of the primary answer, not hedges.
Judge the wording only."""

client = anthropic.Anthropic(timeout=60, max_retries=6)
# primary_answer comes first so the judge writes down the committed cell type (state words included,
# hedges dropped) before choosing a label. Only `label` is used for scoring.
SCHEMA = {"type": "object",
          "properties": {"primary_answer": {"type": "string"}, "label": {"type": "string", "enum": OPTIONS}},
          "required": ["primary_answer", "label"], "additionalProperties": False}


def judge(answer):
    resp = client.messages.create(
        model=JUDGE_MODEL,
        max_tokens=1000,
        system=INSTRUCTIONS,
        messages=[{"role": "user", "content": f"Answer to map: {answer}"}],
        output_config={"effort": "low", "format": {"type": "json_schema", "schema": SCHEMA}},
    )
    text = next(b.text for b in resp.content if b.type == "text")
    out = json.loads(text)
    # Label selection sometimes slips between similar names (primary_answer "NK Proliferating",
    # label "NK_CD56bright"). If the judge's own primary answer is exactly a label or synonym, trust it.
    label = rule_match(out["primary_answer"]) or out["label"]
    # Guard against near-misses like "NK Cell (subtype unclear)": match case-insensitively.
    by_lower = {o.lower(): o for o in OPTIONS}
    if label.lower() not in by_lower:
        raise ValueError(f"judge returned unknown label {label!r} for {answer!r}")
    return by_lower[label.lower()]


def results_files():
    """Model and agent results files only (not scores_hao.csv or the hand-check sheet)."""
    return [p for p in glob.glob("results/*_hao.csv") + glob.glob("results/agents/*_hao.csv")
            if "knob" in pd.read_csv(p, nrows=0).columns]


if __name__ == "__main__":
    fresh = "--fresh" in sys.argv
    cache = {} if fresh or not os.path.exists(CACHE) else json.load(open(CACHE))
    answers = set()
    for path in results_files():
        answers |= set(pd.read_csv(path).predicted.dropna().astype(str))
    ruled = {a for a in answers if rule_match(a)}
    todo = sorted(answers - ruled - set(cache))
    print(f"{len(answers)} distinct answers: {len(ruled)} matched by rules, {len(todo)} to send to the judge")
    with ThreadPoolExecutor(8) as pool:
        for answer, label in zip(todo, pool.map(judge, todo)):
            cache[answer] = label
    json.dump(cache, open(CACHE, "w"), indent=1, sort_keys=True)
    print(f"Saved {len(cache)} judgements to {CACHE}")
