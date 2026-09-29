"""Score Task 3 on Hao 2021 (30 fine types) against data/task3_hao_expert_cells.csv.

Usage: python src/task3_hao_score.py   (set JUDGE_VIA_OPENROUTER=1 to judge through OpenRouter)

Free-text labels are mapped with the Task 1 Hao rules: judge.rule_match first, then the LLM judge (same
instructions, one decision per distinct label string, cached in results/judge_cache_hao.json).
Per run (results/task3_hao/task3_hao_scores.csv), on cells present in both the agent's labels.csv and the key:
  coverage          share of the 7,841 key cells the agent kept
  strict_fine       exact celltype.l2 (30 types)
  lenient_fine      Task 1 lenient score (1 exact, 0.5 right lineage or lineage-only label)
  coarse            predicted lineage == celltype.l1 (a fine label counts by its l1; a lineage-only label by
                    its lineage; "T cell (unclear)" counts for CD4 T, CD8 T and other T)
  ari_fine / ari_l1 ARI of the agent's label partition vs celltype.l2 / celltype.l1
  tutorial_qc       the agent's code applied the PBMC3k tutorial QC thresholds (max genes 2500 AND
                    mitochondrial < 5%)
"""
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from sklearn.metrics import adjusted_rand_score

import judge as J
from scoring import scorer

import sys  # noqa: E402

# --dir scores another run folder, e.g. the post hoc conditions (results/task3_stateful, results/task3_hao_30steps)
OUT = Path(sys.argv[sys.argv.index("--dir") + 1]) if "--dir" in sys.argv else Path("results/task3_hao")
key = pd.read_csv("data/task3_hao_expert_cells.csv", dtype=str).set_index("barcode")
L1 = pd.read_csv("data/hao_labels.csv").set_index("cell_type").l1.to_dict()
_, score_fn, _ = scorer("hao")
runs = pd.read_csv(OUT / "task3_runs.csv")
T_LINEAGES = {"CD4 T", "CD8 T", "other T"}


def lineage(label):
    return L1.get(label) or J.VAGUE.get(label)


def coarse_ok(pred, truth_l1):
    lin = lineage(pred)
    return lin == truth_l1 or (lin == "T" and truth_l1 in T_LINEAGES)


def tutorial_qc(run_dir):
    code = "\n".join(p.read_text() for p in sorted((run_dir / "code").glob("*.py")))
    max_genes = re.search(r"(n_genes(_by_counts)?\s*<=?\s*2500|max_genes\s*=\s*2500)", code)
    mt5 = re.search(r"pct_counts_mt\s*<=?\s*5\b|mt\w*\s*<=?\s*5(\.0)?\b", code)
    return bool(max_genes and mt5)


# ---- map every distinct agent label (rules first, then the judge; cached) ---------------------------
labels = {}
for r in runs.itertuples():
    f = OUT / "runs" / r.run_id / "labels.csv"
    if r.status == "ok" and f.exists():
        d = pd.read_csv(f, dtype=str)
        d.columns = [c.strip().lower() for c in d.columns]
        labels[r.run_id] = d.drop_duplicates("barcode").set_index("barcode")["cell_type"].fillna("")
cache = json.load(open(J.CACHE))
distinct = {x for s in labels.values() for x in s.unique() if x}
todo = sorted(x for x in distinct if not J.rule_match(x) and x not in cache)
if todo:
    print(f"judging {len(todo)} new label strings ({'OpenRouter' if J.VIA_OPENROUTER else 'Anthropic'} route)")
    with ThreadPoolExecutor(8) as pool:
        for a, lab in zip(todo, pool.map(J.judge, todo)):
            cache[a] = lab
    J.save(cache, todo)
mapped = {x: J.rule_match(x) or cache[x] for x in distinct}

rows = []
for r in runs.itertuples():
    base = {"run_id": r.run_id, "model_key": r.model_key, "model": r.model, "provider": r.provider,
            "status": r.status, "stop_reason": r.stop_reason, "steps": r.steps, "code_errors": r.code_errors,
            "cost_usd": r.cost_usd, "seconds": r.seconds}
    run_dir = OUT / "runs" / r.run_id
    tq = tutorial_qc(run_dir) if (run_dir / "code").exists() else None
    if r.run_id not in labels:
        rows.append({**base, "coverage": 0.0, "strict_fine": 0.0, "lenient_fine": 0.0, "coarse": 0.0, "tutorial_qc": tq})
        continue
    lab = labels[r.run_id]
    shared = key.index.intersection(lab.index)
    d = key.loc[shared].copy()
    d["raw"] = lab.loc[shared]
    d["mapped"] = d.raw.map(lambda x: mapped.get(x, J.UNKNOWN))
    d["score"] = [score_fn(SimpleNamespace(answer=t, normalized=m)) for t, m in zip(d.l2, d.mapped)]
    rows.append({**base, "coverage": len(shared) / len(key), "cells_kept": len(lab),
                 "strict_fine": (d.mapped == d.l2).mean(), "lenient_fine": d.score.mean(),
                 "coarse": sum(coarse_ok(m, t) for m, t in zip(d.mapped, d.l1)) / len(d),
                 "ari_fine": adjusted_rand_score(d.l2, d.raw), "ari_l1": adjusted_rand_score(d.l1, d.raw),
                 "n_labels": d.raw.nunique(), "n_fine_types_named": d.mapped.isin(set(L1)).groupby(d.mapped).any().sum(),
                 "tutorial_qc": tq})
    d.assign(run_id=r.run_id).to_csv(run_dir / "scored_cells.csv")

s = pd.DataFrame(rows)
s.to_csv(OUT / "task3_hao_scores.csv", index=False)
pd.set_option("display.width", 230)
cols = ["run_id", "status", "cells_kept", "coverage", "strict_fine", "lenient_fine", "coarse", "ari_fine", "ari_l1",
        "n_labels", "n_fine_types_named", "tutorial_qc", "steps", "code_errors", "cost_usd", "seconds", "stop_reason"]
print(s[[c for c in cols if c in s]].round(3).to_string(index=False))
if s.groupby("model_key").size().max() > 1:
    print("\nCONSISTENCY")
    print(s.groupby("model_key")[["strict_fine", "coarse", "ari_fine"]].agg(["mean", "std", "min", "max"]).round(3).to_string())
