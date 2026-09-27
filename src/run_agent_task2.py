"""Task 2 specialist agent: Claude Sonnet 5 as a LangGraph prebuilt ReAct agent with three tools.

Usage: python src/run_agent_task2.py [--limit N] [--test]

Tools (no tool can see perturbation results - only control cells, MyGene.info and CollecTRI):
  gene_info                    MyGene.info: symbol or Ensembl ID -> symbol, name, summary
  coexpression_in_control      Pearson r of gene X and gene Y across the 11,855 CONTROL K562 cells of
                               Norman 2019 (the same numbers as the co-expression baseline)
  collectri_lookup             CollecTRI TF -> target links (same snapshot as the CollecTRI baseline)

Same setup as the Task 1 agents: claude-sonnet-5 via the Anthropic API, the question text as the only
prompt (no system prompt), no temperature, default thinking, max_tokens 4000, at most 5 tool calls
(later calls get a "budget used up" reply), 2-minute timeout per question.

Outputs: results/agents_task2/specialist/<id>.json (full trajectory) and
         results/agents_task2/agent-specialist_task2.csv (one row per question). Resumable.
--test asks one question, prints the full trajectory, and saves nothing.
"""
import argparse
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout

import numpy as np
import pandas as pd
import requests
from dotenv import load_dotenv
from langchain_anthropic import ChatAnthropic
from langchain_core.messages import AIMessage, HumanMessage, ToolMessage
from langchain_core.tools import StructuredTool
from langgraph.errors import GraphRecursionError
from langgraph.prebuilt import create_react_agent

from prices import cost_usd
from task2_common import parse_direction

load_dotenv()
parser = argparse.ArgumentParser()
parser.add_argument("--limit", type=int, help="pilot: only N questions, spread evenly over the set")
parser.add_argument("--test", action="store_true", help="one question, print the trajectory, save nothing")
args = parser.parse_args()

MODEL, MAX_TOOL_CALLS, TIMEOUT_S, WORKERS = "claude-sonnet-5", 5, 120, 4
OUT_DIR, CSV = "results/agents_task2/specialist", "results/agents_task2/agent-specialist_task2.csv"

# ---- tool data ------------------------------------------------------------------------------------
_c = np.load("data/norman2019/task2_ctrl_corr.npz", allow_pickle=True)
CORR, X_GENES, GENES = _c["corr"], list(_c["x_genes"]), list(_c["genes"])
X_IDX, G_IDX = {g: i for i, g in enumerate(X_GENES)}, {g: i for i, g in enumerate(GENES)}
# Symbols are case-sensitive (e.g. C1orf56); match exactly first, then case-insensitively.
_CANON = {g.upper(): g for g in GENES}
_CANON.update({g.upper(): g for g in X_GENES})


def canon(symbol):
    s = symbol.strip()
    return s if s in G_IDX or s in X_IDX else _CANON.get(s.upper(), s)
CTRL_FRAC, X_CTRL_FRAC = _c["ctrl_frac"], _c["x_ctrl_frac"]
NET = pd.read_csv("data/task2/collectri_human.csv")


def gene_info(genes: list[str]) -> str:
    """Look up genes on MyGene.info. Accepts gene symbols or Ensembl gene IDs (ENSG...), up to 20 at
    a time. Returns each gene's official symbol, full name and a short functional summary."""
    r = requests.post("https://mygene.info/v3/query", timeout=30, data={
        "q": ",".join(genes[:20]), "scopes": "symbol,ensembl.gene", "species": "human",
        "fields": "symbol,name,summary"})
    r.raise_for_status()
    out = {}
    for hit in r.json():  # several hits per query are possible; prefer one with a summary
        q = hit["query"]
        if q in out and "notfound" not in out[q] and (out[q]["summary"] or not hit.get("summary")):
            continue
        out[q] = ({"notfound": True} if hit.get("notfound") else
                  {"symbol": hit.get("symbol"), "name": hit.get("name"), "summary": (hit.get("summary") or "")[:400]})
    return json.dumps(out, separators=(",", ":"))


def coexpression_in_control(gene_x: str, gene_y: str) -> str:
    """Pearson correlation of two genes' expression across unperturbed (control) K562 cells from the
    Norman 2019 CRISPRa screen, with the share of control cells expressing each gene. gene_x must be one
    of the genes activated in this screen; gene_y can be any measured gene. Uses control cells only."""
    x, y = canon(gene_x), canon(gene_y)
    if x not in X_IDX and y in X_IDX:  # accept the pair in either order
        x, y = y, x
    if x not in X_IDX:
        return json.dumps({"error": f"{gene_x} is not one of the activated genes in this screen"})
    if y not in G_IDX:
        return json.dumps({"error": f"{gene_y} is not measured"})
    r = float(CORR[X_IDX[x], G_IDX[y]])
    return json.dumps({"gene_x": x, "gene_y": y, "pearson_r_control_cells": None if np.isnan(r) else round(r, 4),
                       "n_control_cells": 11855, "share_control_cells_expressing_x": round(float(X_CTRL_FRAC[X_IDX[x]]), 4),
                       "share_control_cells_expressing_y": round(float(CTRL_FRAC[G_IDX[y]]), 4)})


def collectri_lookup(tf: str, target: str = "") -> str:
    """Look up the CollecTRI transcription-factor -> target network (curated TF regulons). Returns
    whether `tf` is in CollecTRI, whether it has a link to `target` (with sign: +1 activates, -1
    represses), and up to 30 of its other targets. Links without sign evidence are marked
    "default activation"."""
    tf, target = canon(tf), canon(target) if target else ""
    sub = NET[NET.source == tf]
    out = {"tf": tf, "tf_in_collectri": bool(len(sub)), "n_targets": int(len(sub))}
    if target:
        hit = sub[sub.target == target]
        out["link_to_target"] = None if hit.empty else {
            "target": target, "weight": int(hit.weight.iloc[0]), "sign_decision": hit.sign_decision.iloc[0],
            "n_references": len(str(hit.references.iloc[0]).split(";"))}
    out["some_targets"] = [{"target": t, "weight": int(w)} for t, w in zip(sub.target.head(30), sub.weight.head(30))]
    return json.dumps(out, separators=(",", ":"))


TOOLS = [gene_info, coexpression_in_control, collectri_lookup]


def make_tools(counter):
    """Fresh tools per question sharing one call counter, so the budget covers all tools together."""
    def wrap(fn):
        def run(**kwargs):
            with counter["lock"]:
                counter["n"] += 1
                n = counter["n"]
            if n > MAX_TOOL_CALLS:
                return (f"Tool budget used up ({MAX_TOOL_CALLS} calls). Do not call tools again; "
                        "give your final JSON answer now.")
            try:
                return fn(**kwargs)
            except Exception as e:
                return f"Tool error: {type(e).__name__}: {str(e)[:200]}"
        return StructuredTool.from_function(func=run, name=fn.__name__, description=fn.__doc__,
                                            args_schema=StructuredTool.from_function(fn).args_schema)
    return [wrap(f) for f in TOOLS]


llm = ChatAnthropic(model=MODEL, max_tokens=4000, timeout=TIMEOUT_S, max_retries=6)


def serialize(m):
    d = {"type": m.type, "content": m.content}
    if isinstance(m, AIMessage):
        d["tool_calls"] = [{"name": c["name"], "args": c["args"], "id": c["id"]} for c in m.tool_calls]
        d["usage"], d["stop_reason"] = m.usage_metadata, (m.response_metadata or {}).get("stop_reason")
    if isinstance(m, ToolMessage):
        d["tool_name"], d["tool_call_id"] = m.name, m.tool_call_id
    return d


def final_text(messages):
    for m in reversed(messages):
        if isinstance(m, AIMessage):
            c = m.content
            return c if isinstance(c, str) else "".join(b.get("text", "") for b in c if isinstance(b, dict))
    return ""


def ask(q):
    counter = {"n": 0, "lock": threading.Lock()}
    agent = create_react_agent(llm, make_tools(counter))
    result, status, start = {"messages": [HumanMessage(q["question"])]}, "ok", time.monotonic()
    pool = ThreadPoolExecutor(1)
    fut = pool.submit(agent.invoke, {"messages": [HumanMessage(q["question"])]},
                      {"recursion_limit": 2 * MAX_TOOL_CALLS + 6})
    try:
        result = fut.result(timeout=TIMEOUT_S)
    except FuturesTimeout:
        status = "timeout"
    except GraphRecursionError:
        status = "recursion_limit"
    except Exception as e:  # e.g. network down after retries: not saved, so a rerun retries it
        status = f"error: {type(e).__name__}: {str(e)[:120]}"
    pool.shutdown(wait=False)
    messages = result["messages"]
    raw = final_text(messages) if status == "ok" else ""
    pred, conf = parse_direction(raw)
    ai = [m for m in messages if isinstance(m, AIMessage)]
    tin = sum((m.usage_metadata or {}).get("input_tokens", 0) for m in ai)
    tout = sum((m.usage_metadata or {}).get("output_tokens", 0) for m in ai)
    row = {"id": q["id"], "X": q["X"], "Y": q["Y"], "label": q["label"], "predicted": pred, "confidence": conf,
           "raw": raw, "finish_reason": status, "prompt_tokens": tin, "completion_tokens": tout,
           "cost_usd": cost_usd(MODEL, tin, tout), "n_tool_calls": counter["n"], "n_model_calls": len(ai),
           "latency_s": round(time.monotonic() - start, 1)}
    traj = {"condition": "task2_specialist", "model": MODEL, "question": q, "status": status,
            "messages": [serialize(m) for m in messages], "summary": row}
    return row, traj


questions = json.load(open("data/task2_questions.json"))
if args.test:
    row, traj = ask(questions[0])
    print(json.dumps(traj, indent=1, default=str)[:6000])
    raise SystemExit
os.makedirs(OUT_DIR, exist_ok=True)
todo_q = questions[:: len(questions) // args.limit][: args.limit] if args.limit else questions
done = pd.read_csv(CSV) if os.path.exists(CSV) else pd.DataFrame()
todo = [q for q in todo_q if done.empty or q["id"] not in set(done.id)]
print(f"task2 specialist: {len(todo)} to ask, {len(done)} already done", flush=True)
rows, failed = [], 0
with ThreadPoolExecutor(WORKERS) as pool:
    for i, (row, traj) in enumerate(pool.map(ask, todo), 1):
        if row["finish_reason"].startswith("error"):
            failed += 1
            print(f"[{i}/{len(todo)}] {row['id']} FAILED (will retry on rerun): {row['finish_reason']}", flush=True)
            continue
        json.dump(traj, open(f"{OUT_DIR}/{row['id']}.json", "w"), indent=1, default=str)
        rows.append(row)
        print(f"[{i}/{len(todo)}] {row['id']} {row['X']}->{row['Y']} truth={row['label']:<9} pred={row['predicted']} "
              f"tools={row['n_tool_calls']} ${row['cost_usd']:.4f}", flush=True)
        results = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
        results.to_csv(CSV + ".tmp", index=False)
        os.replace(CSV + ".tmp", CSV)
print(f"\nSaved {len(done) + len(rows)} rows to {CSV} | this run ${sum(r['cost_usd'] for r in rows):.3f} | failed: {failed}")
