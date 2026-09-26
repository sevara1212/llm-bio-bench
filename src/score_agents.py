"""Score the agent conditions against plain Claude and the CellMarker lookup, on the same 96 questions.

Run judge.py first (it judges agent answers too). Uses scoring.py - the same label mapping and scores
as everything else.

Prints accuracy per condition (strict, lenient) overall, by gene format and by lineage group, then
classifies every wrong specialist answer:
  no_tool            answered without calling any tool
  no_cellmarker      called only gene_info (e.g. converted Ensembl IDs), never queried CellMarker
  cellmarker_lacked  queried CellMarker, but no output mentioned a cell type mapping to the truth
  truth_tied_top     the true type was among several types sharing the most support in CellMarker output
  truth_below_top    the true type appeared in CellMarker output, but another type had more support
  ignored            the true type alone had the most support in CellMarker output, yet the agent answered otherwise
  no_answer          timeout / no parsable answer
Support is counted from cellmarker_gene_to_cell_types outputs only.
"Support" = number of queried genes whose CellMarker output lists a type mapping to that label (the
same counting idea as the lookup baseline). For wrong answers the scorer also says if it got lineage
credit (lenient 0.5) and, for Ensembl questions, whether gene_info returned the right symbols.

Writes results/agents/scores_agents_hao.csv and results/agents/specialist_errors_hao.csv.
"""
import csv
import json
from collections import Counter
from types import SimpleNamespace

import pandas as pd

from scoring import scorer

pd.set_option("display.width", 220)
pd.set_option("display.max_colwidth", 70)
normalize, score, _ = scorer("hao")
ids = json.load(open("data/agent_subset_hao.json"))["ids"]
L1 = pd.read_csv("data/hao_labels.csv").set_index("cell_type").l1
CM_TO_HAO = dict(csv.reader(open("configs/cellmarker_to_hao.csv")))
ID_TO_SYMBOL = {e: g for f in ["data/hao_markers.csv", "data/hao_background_genes.csv"]
                for e, g in zip(pd.read_csv(f).ensembl, pd.read_csv(f).gene)}

SOURCES = {
    "plain Claude": "results/claude-sonnet-5_hao.csv",
    "CellMarker lookup": "results/cellmarker-lookup_hao.csv",
    "generic agent": "results/agents/agent-generic_hao.csv",
    "specialist agent": "results/agents/agent-specialist_hao.csv",
}
frames = []
for name, path in SOURCES.items():
    d = pd.read_csv(path)
    d = d[d.id.isin(ids) & (d.run == 0)].copy()
    missing = set(ids) - set(d.id)
    assert not missing, f"{name}: {len(missing)} subset questions missing, e.g. {sorted(missing)[:3]}"
    d.insert(0, "condition", name)
    frames.append(d)
s = pd.concat(frames, ignore_index=True)
s["normalized"] = s.predicted.map(normalize)
not_judged = s.normalized == "NOT JUDGED"
assert not not_judged.any(), f"{not_judged.sum()} answers not judged yet: run python src/judge.py"
s["score"] = [score(SimpleNamespace(answer=a, normalized=n)) for a, n in zip(s.answer, s.normalized)]
s["strict"] = (s.score == 1).astype(float)
s["lineage"] = s.answer.map(L1)
s.to_csv("results/agents/scores_agents_hao.csv", index=False)

order = list(SOURCES)
n = s.groupby("condition").size()
print(f"n = {n.iloc[0]} questions per condition ({(s[s.condition == order[0]].format == 'symbol').sum()} symbol, "
      f"{(s[s.condition == order[0]].format == 'ensembl').sum()} Ensembl)\n")
t = s.groupby("condition").agg(strict=("strict", "mean"), lenient=("score", "mean")).reindex(order)
print("OVERALL"); print(t.round(3).to_string())
t = s.pivot_table(index="condition", columns="format", values=["strict", "score"]).reindex(order)
t.columns = [f"{'lenient' if m == 'score' else m} ({f})" for m, f in t.columns]
print("\nBY GENE FORMAT (n = 48 each)"); print(t.round(3).to_string())
t = s.pivot_table(index="lineage", columns="condition", values="strict")[order]
t.insert(0, "n", s[s.condition == order[0]].groupby("lineage").size())
print("\nSTRICT BY LINEAGE GROUP"); print(t.round(2).to_string())
t = s.pivot_table(index="lineage", columns="condition", values="score")[order]
print("\nLENIENT BY LINEAGE GROUP"); print(t.round(2).to_string())

print("\nPAIRED, SAME 96 QUESTIONS (strict): specialist vs each other condition")
piv = s.pivot_table(index="id", columns="condition", values="strict")
for other in ["plain Claude", "CellMarker lookup", "generic agent"]:
    gained = int(((piv["specialist agent"] == 1) & (piv[other] == 0)).sum())
    lost = int(((piv["specialist agent"] == 0) & (piv[other] == 1)).sum())
    print(f"   vs {other:<18} specialist right & other wrong: {gained:>2} | other right & specialist wrong: {lost:>2}")

agents = s[s.condition.str.contains("agent")]
print("\nAGENT BEHAVIOUR")
print(agents.groupby("condition").agg(
    mean_tool_calls=("n_tool_calls", "mean"), share_no_tool=("n_tool_calls", lambda x: (x == 0).mean()),
    share_hit_budget=("n_tool_calls", lambda x: (x >= 5).mean()), not_ok=("finish_reason", lambda x: (x != "ok").sum()),
    mean_latency_s=("latency_s", "mean"), total_cost=("cost_usd", "sum"), cost_per_q=("cost_usd", "mean")).round(3).to_string())


# ---- classify the specialist's wrong answers ----------------------------------------------------
def tool_evidence(traj):
    """(tools used, label support Counter from CellMarker outputs, gene_info symbols ok?)."""
    used, support, symbols_ok = [], Counter(), None
    calls = {c["id"]: c for m in traj["messages"] if m["type"] == "ai" for c in m.get("tool_calls", [])}
    for m in traj["messages"]:
        if m["type"] != "tool":
            continue
        call = calls.get(m["tool_call_id"], {})
        used.append(m["tool_name"])
        try:
            out = json.loads(m["content"])
        except (json.JSONDecodeError, TypeError):
            continue  # budget message or tool error
        if m["tool_name"] == "cellmarker_gene_to_cell_types":
            for gene, info in out.items():
                labels = {CM_TO_HAO.get(x["cell_type"]) for x in info["cell_types"]} - {None}
                support.update(labels)
        elif m["tool_name"] == "cellmarker_cell_type_markers" and "markers" in out:
            asked = {str(ID_TO_SYMBOL.get(g, g)).upper() for g in traj["question"]["genes"]}
            if asked & set(out["markers"]):
                support[CM_TO_HAO.get(out["cell_type"])] += 0  # looked at, but not counted as gene support
        elif m["tool_name"] == "gene_info":
            # Only IDs the dataset itself names (some Ensembl IDs have no symbol in the dataset).
            wanted = [g for g in call.get("args", {}).get("genes", [])
                      if g in ID_TO_SYMBOL and not str(ID_TO_SYMBOL[g]).startswith("ENSG")]
            wrong = [f"{g}->{(out.get(g) or {}).get('symbol')} (dataset: {ID_TO_SYMBOL[g]})"
                     for g in wanted if (out.get(g) or {}).get("symbol") != ID_TO_SYMBOL[g]]
            if wanted:
                symbols_ok = "; ".join(wrong) if wrong else True
    return used, support, symbols_ok


spec = s[s.condition == "specialist agent"]
rows = []
for r in spec[spec.strict < 1].itertuples():
    traj = json.load(open(f"results/agents/specialist/{r.id}.json"))
    used, support, symbols_ok = tool_evidence(traj)
    truth_support = support.get(r.answer, 0)
    top = max(support.values()) if support else 0
    if r.finish_reason != "ok" or not isinstance(r.predicted, str):
        cls = "no_answer"
    elif not used:
        cls = "no_tool"
    elif "cellmarker_gene_to_cell_types" not in used:
        cls = "no_cellmarker"
    elif truth_support == 0:
        cls = "cellmarker_lacked"
    elif truth_support == top and sum(v == top for v in support.values()) == 1:
        cls = "ignored"
    elif truth_support == top:
        cls = "truth_tied_top"
    else:
        cls = "truth_below_top"
    rows.append({"id": r.id, "format": r.format, "truth": r.answer, "predicted": r.predicted, "mapped_to": r.normalized,
                 "lenient": r.score, "class": cls, "tools_used": ",".join(used) or "-",
                 "truth_support": truth_support, "top_support": top,
                 "top_labels": ", ".join(k for k, v in support.most_common(3)),
                 "gene_info_symbols_ok": symbols_ok})
errs = pd.DataFrame(rows)
errs.to_csv("results/agents/specialist_errors_hao.csv", index=False)
print(f"\nSPECIALIST WRONG ANSWERS (strict): {len(errs)} of {len(spec)}")
print(errs["class"].value_counts().to_string())
print("\nwith lineage credit (lenient 0.5):", int((errs.lenient == 0.5).sum()))
wrong_sym = errs[errs.gene_info_symbols_ok.apply(lambda x: isinstance(x, str))]
print(f"gene_info returned a different symbol than the dataset on {len(wrong_sym)} wrong answers:")
for r in wrong_sym.itertuples():
    print(f"   {r.id}: {r.gene_info_symbols_ok}")
print("\nby format:"); print(errs.groupby(["class", "format"]).size().unstack(fill_value=0).to_string())
print()
print(errs[["id", "truth", "predicted", "class", "tools_used", "truth_support", "top_support", "top_labels"]].to_string(index=False))
