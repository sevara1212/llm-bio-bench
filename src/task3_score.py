"""Score Task 3 runs against the expert PBMC3k labels (data/pbmc3k_expert_cells.csv).

Usage: python src/task3_score.py

Per run (results/task3/task3_scores.csv):
  coverage        share of the 2,638 expert cells present in the agent's labels.csv; only shared cells are
                  scored. extra_cells = agent cells the expert pipeline removed in QC.
  strict/lenient  per-cell accuracy after mapping the agent's free-text labels to the 8 PBMC3k types with
                  the Task 1 PBMC3k rules (scoring.py: normalize_pbmc / score_pbmc)
  ARI             adjusted Rand index of the agent's label partition vs the expert Leiden clusters; also
                  ARI of the agent's own cluster column when it saved one in an .h5ad (ari_agent_clusters)
  n_labels        distinct cell_type labels; n_agent_clusters from the agent's saved clusters if any
  steps, code_errors, cost_usd, seconds from task3_runs.csv
Automatic first-pass diagnosis per expert cell type scored wrong (results/task3/task3_diagnosis.csv):
  qc_strict / qc_loose       coverage < 90% / extra cells > 5% of the expert set
  merged_in_clustering       most of the type's cells sit in agent clusters whose majority is another type
  separate_cluster_wrong_name the type has its own agent cluster(s) but the label maps to something else
  no_cluster_info            the agent saved no clusters, so the step can't be read off automatically
Consistency per model: mean, SD, min and max over runs of strict accuracy, ARI and n_labels.
"""
import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

from scoring import normalize_pbmc, score_pbmc

import sys  # noqa: E402

# --dir scores another run folder, e.g. the post hoc conditions (results/task3_stateful, results/task3_hao_30steps)
OUT = Path(sys.argv[sys.argv.index("--dir") + 1]) if "--dir" in sys.argv else Path("results/task3")
expert = pd.read_csv("data/pbmc3k_expert_cells.csv", dtype=str).set_index("barcode")
runs = pd.read_csv(OUT / "task3_runs.csv")
CLUSTER_COLS = ("leiden", "louvain", "cluster", "clusters")


def agent_clusters(run_dir, labels):
    """The agent's own cluster assignment: among the cluster-like columns of every .h5ad obs table it
    saved, the one that best matches its labels.csv (highest ARI) - i.e. the clustering it labelled from,
    not e.g. a resolution it tried and discarded."""
    best, best_ari = None, -1.0
    for f in sorted(run_dir.glob("agent_obs_*.csv")):
        obs = pd.read_csv(f, index_col=0, dtype=str)
        for col in obs.columns:
            if any(k in col.lower() for k in CLUSTER_COLS) and 2 <= obs[col].nunique() <= 60:
                common = obs.index.intersection(labels.index)
                if len(common) < 100:
                    continue
                ari = adjusted_rand_score(labels.loc[common], obs.loc[common, col])
                if ari > best_ari:
                    best, best_ari = obs[col].rename(f"{f.stem}:{col}"), ari
    return best


rows, diag = [], []
for r in runs.itertuples():
    run_dir = OUT / "runs" / r.run_id
    base = {"run_id": r.run_id, "model_key": r.model_key, "model": r.model, "provider": r.provider,
            "status": r.status, "stop_reason": r.stop_reason, "steps": r.steps, "code_errors": r.code_errors,
            "cost_usd": r.cost_usd, "seconds": r.seconds}
    f = run_dir / "labels.csv"
    if r.status != "ok" or not f.exists():
        rows.append({**base, "coverage": 0.0, "strict": 0.0, "lenient": 0.0})
        diag.append({"run_id": r.run_id, "issue": "ran_out_or_failed", "detail": r.stop_reason})
        continue
    lab = pd.read_csv(f, dtype=str)
    lab.columns = [c.strip().lower() for c in lab.columns]
    if not {"barcode", "cell_type"} <= set(lab.columns):
        rows.append({**base, "status": "bad_labels_file", "coverage": 0.0, "strict": 0.0, "lenient": 0.0})
        diag.append({"run_id": r.run_id, "issue": "bad_labels_file", "detail": f"columns {list(lab.columns)}"})
        continue
    lab = lab.drop_duplicates("barcode").set_index("barcode")
    shared = expert.index.intersection(lab.index)
    d = expert.loc[shared].copy()
    d["raw"] = lab.loc[shared, "cell_type"]
    d["normalized"] = d.raw.map(normalize_pbmc)
    d["score"] = [score_pbmc(pd.Series({"normalized": n, "answer": a})) for n, a in zip(d.normalized, d.cell_type)]
    clusters = agent_clusters(run_dir, lab.cell_type)
    ari_cl = None
    if clusters is not None:
        common = shared.intersection(clusters.index)
        ari_cl = adjusted_rand_score(d.loc[common, "cluster"], clusters.loc[common]) if len(common) else None
    coverage, extra = len(shared) / len(expert), len(lab.index.difference(expert.index))
    rows.append({**base, "coverage": coverage, "cells_kept": len(lab), "extra_cells": extra,
                 "strict": (d.score == 1).mean(), "lenient": d.score.mean(),
                 "ari_labels": adjusted_rand_score(d.cluster, d.raw), "ari_agent_clusters": ari_cl,
                 "n_labels": d.raw.nunique(), "n_agent_clusters": None if clusters is None else clusters.nunique(),
                 "cluster_column": None if clusters is None else clusters.name})
    # ---- first-pass diagnosis --------------------------------------------------------------------
    if coverage < 0.90:
        diag.append({"run_id": r.run_id, "issue": "qc_strict", "detail": f"coverage {coverage:.1%}"})
    if extra > 0.05 * len(expert):
        diag.append({"run_id": r.run_id, "issue": "qc_loose", "detail": f"{extra} cells the expert QC removed"})
    if clusters is not None:  # majority expert type of each agent cluster
        cl_all = clusters.reindex(d.index).dropna()
        majority = d.loc[cl_all.index].groupby(cl_all).cell_type.agg(lambda x: x.value_counts().index[0])
    for ct, g in d.groupby("cell_type"):
        wrong = g[g.score < 1]
        if len(wrong) < 0.2 * len(g):  # >= 80% right: not flagged
            continue
        wl = wrong.raw.value_counts()
        detail = (f"{ct}: {1 - len(wrong) / len(g):.0%} strict; {len(wrong)}/{len(g)} cells wrong, mostly labelled "
                  f"{wl.index[0]!r} -> {normalize_pbmc(wl.index[0])}")
        if clusters is None:
            issue = "no_cluster_info"
        else:
            wc = cl_all.reindex(wrong.index).dropna()
            own = wc.map(majority).eq(ct).mean() if len(wc) else 0.0
            # wrong cells in clusters where this type is the majority -> the cluster was named wrong;
            # wrong cells in clusters dominated by another type -> they were merged into that cluster
            issue = "separate_cluster_wrong_name" if own >= 0.5 else "merged_in_clustering"
            if issue == "separate_cluster_wrong_name":
                # the type is the majority of its cluster - but is that cluster really pure, or does it also
                # hold >= 25% of another type (a merged cluster named after the minority type)?
                top_cl = wc[wc.map(majority).eq(ct)].value_counts().index[0]
                mix = d.loc[cl_all.index[cl_all == top_cl], "cell_type"].value_counts(normalize=True)
                if len(mix) > 1 and mix.iloc[1] >= 0.25:
                    issue = "merged_in_clustering"
                    detail += (f"; its cluster {top_cl} is a mix ({mix.index[0]} {mix.iloc[0]:.0%}, "
                               f"{mix.index[1]} {mix.iloc[1]:.0%}) named after the minority type")
            other = wc.map(majority).value_counts()
            detail += f"; {own:.0%} of the wrong cells sit in clusters where {ct} is the majority"
            if issue == "merged_in_clustering" and len(other):
                detail += f" (merged into a {other.index[0]}-majority cluster)"
        diag.append({"run_id": r.run_id, "issue": issue, "detail": detail})

scores = pd.DataFrame(rows)
scores.to_csv(OUT / "task3_scores.csv", index=False)
pd.DataFrame(diag).to_csv(OUT / "task3_diagnosis.csv", index=False)
pd.set_option("display.width", 220)
cols = ["run_id", "status", "coverage", "strict", "lenient", "ari_labels", "ari_agent_clusters", "n_labels",
        "n_agent_clusters", "cells_kept", "steps", "code_errors", "cost_usd", "seconds", "stop_reason"]
print(scores[[c for c in cols if c in scores]].round(3).to_string(index=False))
if len(scores.groupby("model_key")) and scores.groupby("model_key").size().max() > 1:
    print("\nCONSISTENCY (per model, over runs)")
    print(scores.groupby("model_key")[["strict", "ari_labels", "n_labels"]].agg(["mean", "std", "min", "max"]).round(3).to_string())
print("\nFIRST-PASS DIAGNOSIS")
print(pd.DataFrame(diag).to_string(index=False) if diag else "(no issues found)")
