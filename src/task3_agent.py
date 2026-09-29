"""Task 3: an agent analyses raw PBMC3k counts end to end in a sandbox (QC -> labels.csv).

Usage: python src/task3_agent.py --runs N [--models claude,gpt,gemini] [--dataset pbmc3k|hao] [--min-balance 3]
  Brings every model up to N runs (runs already done are kept), one run at a time, saving after each.

Per run:
  - fresh sandbox folder with only the raw 10x counts (task3_sandbox.py)
  - LangGraph prebuilt ReAct agent, one tool run_python(code); model via OpenRouter (all three models)
  - the same user instruction for every model (INSTRUCTION); no system prompt, no temperature set
  - at most 15 code executions; stop at $1.50 per run (checked after every model reply)
  - no labels.csv at the end -> failed
Budget: before each run the live OpenRouter balance is checked; everything stops if it is below $3 (or
would fall below $3 during the run).
Saves results/task3/runs/<run_id>/: trajectory.json (every code cell, its output, the model's text and
reasoning where returned, tokens, cost), labels.csv, the executed code cells, and the obs table of any
.h5ad the agent wrote (its own clusters, used for the failure analysis). One row per run is appended to
results/task3/task3_runs.csv.
"""
import argparse
import json
import os
import re
import shutil
import time
import urllib.request
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from langchain_core.messages import AIMessage, HumanMessage, ToolMessage
from langchain_core.tools import StructuredTool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent

from task3_sandbox import StatefulSession, new_run_folder, run_python

load_dotenv()
MODELS = {"claude": "anthropic/claude-sonnet-5", "gpt": "openai/gpt-5.6-terra", "gemini": "google/gemini-3.8-flash"}
PROVIDER = "openrouter"
MAX_STEPS, RUN_BUDGET, MIN_BALANCE, MAX_TOKENS = 15, 1.50, 3.00, 8000
OUT = Path("results/task3")  # set per dataset in main (results/task3 or results/task3_hao)
RUNS_CSV = OUT / "task3_runs.csv"
DATASET = "pbmc3k"
STATEFUL = False  # post hoc condition C: variables persist between run_python calls
STATEFUL_DOC = """Run Python code in a persistent Python session in the working folder (like a Jupyter notebook) and
        return its printed output (stdout and stderr, truncated to about 50 lines / 4,000 characters). Variables,
        imports and loaded data persist between calls. Each call has a 60-second limit; if it times out, the
        session restarts and all variables are lost. scanpy, anndata, pandas, numpy, scipy and leidenalg are
        installed. No internet access."""
INSTRUCTION = ("This folder contains a raw single-cell RNA-seq count matrix from human peripheral blood "
               "mononuclear cells (PBMCs). Perform standard quality control, normalisation and clustering, "
               "identify marker genes, and assign a cell type label to every cluster. Save labels.csv with "
               "columns barcode,cell_type. Work step by step and check your outputs.")
KEY = os.environ["OPENROUTER_API_KEY"]


def openrouter(path):
    req = urllib.request.Request(f"https://openrouter.ai/api/v1/{path}", headers={"Authorization": f"Bearer {KEY}"})
    return json.load(urllib.request.urlopen(req, timeout=30))["data"]


def balance(tries=6):
    """Live OpenRouter balance; retries on network errors (the connection dropped once mid-run)."""
    for i in range(tries):
        try:
            d = openrouter("credits")
            return d["total_credits"] - d["total_usage"]
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(20 * (i + 1))


PRICES = {m["id"]: (float(m["pricing"]["prompt"]), float(m["pricing"]["completion"]))
          for m in openrouter("models") if m["id"] in MODELS.values()}


class Stop(Exception):
    pass


def one_run(model_key, run_id):
    model = MODELS[model_key]
    start_balance = balance()
    if start_balance < MIN_BALANCE:
        raise Stop(f"OpenRouter balance ${start_balance:.2f} is below ${MIN_BALANCE:.2f}")
    run_dir = OUT / "runs" / run_id
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    work = new_run_folder(run_id, DATASET)
    session = StatefulSession(work) if STATEFUL else None
    steps, cells = {"n": 0}, []

    def tool_run_python(code: str) -> str:
        """Run Python code in the working folder and return its printed output (stdout and stderr,
        truncated to about 50 lines / 4,000 characters). Each call runs in a NEW Python process with a
        60-second limit: variables do not persist between calls, so save intermediate results to files in
        the working folder (e.g. .h5ad) and load them in the next call. scanpy, anndata, pandas, numpy,
        scipy and leidenalg are installed. No internet access."""
        if steps["n"] >= MAX_STEPS:
            return f"No code executions left ({MAX_STEPS}/{MAX_STEPS} used). The run is over."
        steps["n"] += 1
        code = re.sub(r"^```(?:python)?\s*|\s*```$", "", code.strip())
        r = session.run(code, steps["n"]) if session else run_python(work, code, steps["n"])
        r["step"] = steps["n"]
        cells.append(r)
        left = MAX_STEPS - steps["n"]
        note = f"[execution {steps['n']}/{MAX_STEPS}, exit code {r['returncode']}; {left} left]"
        if left == 0:
            note += " This was your last code execution."
        return f"{r['output_shown']}\n{note}"

    tool = StructuredTool.from_function(func=tool_run_python, name="run_python",
                                        description=STATEFUL_DOC if STATEFUL else tool_run_python.__doc__)
    llm = ChatOpenAI(model=model, base_url="https://openrouter.ai/api/v1", api_key=KEY, max_tokens=MAX_TOKENS,
                     timeout=180, max_retries=3, extra_body={"usage": {"include": True}})
    agent = create_react_agent(llm, [tool])
    p_in, p_out = PRICES[model]
    messages, cost, tin, tout, stop_reason = [HumanMessage(INSTRUCTION)], 0.0, 0, 0, "agent finished"
    t0 = time.monotonic()
    try:
        for update in agent.stream({"messages": [HumanMessage(INSTRUCTION)]},
                                   {"recursion_limit": 2 * MAX_STEPS + 12}, stream_mode="updates"):
            for node in update.values():
                for m in node.get("messages", []):
                    messages.append(m)
                    if isinstance(m, AIMessage) and m.usage_metadata:
                        tin += m.usage_metadata.get("input_tokens", 0)
                        tout += m.usage_metadata.get("output_tokens", 0)
                        cost = tin * p_in + tout * p_out
            if cost >= RUN_BUDGET:
                stop_reason = f"run budget ${RUN_BUDGET:.2f} reached"
                break
            if start_balance - cost < MIN_BALANCE:
                stop_reason = f"account balance would drop below ${MIN_BALANCE:.0f}"
                break
            if steps["n"] >= MAX_STEPS and isinstance(messages[-1], ToolMessage) and "No code executions left" in str(messages[-1].content):
                stop_reason = "step limit reached"
                break
    except Exception as e:  # recursion limit, API error, ...
        stop_reason = f"error: {type(e).__name__}: {str(e)[:200]}"
    seconds = time.monotonic() - t0
    if session:
        session.close()
    if steps["n"] >= MAX_STEPS and stop_reason == "agent finished":
        stop_reason = "agent finished (all steps used)"

    labels_ok = (work / "labels.csv").exists()
    if labels_ok:
        shutil.copy(work / "labels.csv", run_dir / "labels.csv")
    (run_dir / "code").mkdir()
    for c in work.glob("_cell_*.py"):  # (stateful runs also leave _cell_NN.py.out files; code only is copied)
        shutil.copy(c, run_dir / "code" / c.name)
    for h in work.glob("*.h5ad"):  # the agent's own clusters, for the failure analysis
        try:
            import anndata as ad
            obs = ad.read_h5ad(h, backed="r").obs
            obs.to_csv(run_dir / f"agent_obs_{h.stem}.csv")
        except Exception:
            pass
    files = sorted(str(p.relative_to(work)) for p in work.rglob("*") if p.is_file() and ".tmp" not in p.parts)
    # Free disk space: the agent's intermediate .h5ad files can be GBs per run (Hao). Everything scoring needs
    # (labels.csv, code, obs tables of every .h5ad, the file list) has been copied to results/ above.
    for h in work.rglob("*.h5ad"):
        h.unlink()
    shutil.rmtree(work / ".tmp", ignore_errors=True)

    def ser(m):
        d = {"type": m.type, "content": m.content}
        if isinstance(m, AIMessage):
            d["tool_calls"] = [{"name": c["name"], "args": c["args"]} for c in m.tool_calls]
            d["usage"] = m.usage_metadata
            d["reasoning"] = (m.additional_kwargs or {}).get("reasoning") or (m.additional_kwargs or {}).get("reasoning_content")
        return d
    row = {"run_id": run_id, "dataset": DATASET, "stateful": STATEFUL, "max_steps": MAX_STEPS,
           "model_key": model_key, "model": model, "provider": PROVIDER,
           "status": "ok" if labels_ok else "failed", "stop_reason": stop_reason, "steps": steps["n"],
           "code_errors": sum(c["error"] for c in cells), "timeouts": sum(c["timed_out"] for c in cells),
           "prompt_tokens": tin, "completion_tokens": tout, "cost_usd": round(cost, 4),
           "balance_before": round(start_balance, 2), "seconds": round(seconds, 1)}
    json.dump({"run": row, "instruction": INSTRUCTION, "cells": cells, "files_in_folder": files,
               "messages": [ser(m) for m in messages]}, open(run_dir / "trajectory.json", "w"), indent=1, default=str)
    runs = pd.read_csv(RUNS_CSV) if RUNS_CSV.exists() else pd.DataFrame()
    runs = pd.concat([runs[runs.run_id != run_id] if len(runs) else runs, pd.DataFrame([row])], ignore_index=True)
    runs.to_csv(RUNS_CSV, index=False)
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, required=True)
    ap.add_argument("--models", default="claude,gpt,gemini")
    ap.add_argument("--dataset", default="pbmc3k", choices=["pbmc3k", "hao"])
    ap.add_argument("--min-balance", type=float, default=MIN_BALANCE)
    ap.add_argument("--max-steps", type=int, default=MAX_STEPS, help="post hoc condition B uses 30")
    ap.add_argument("--stateful", action="store_true", help="post hoc condition C: persistent Python session")
    args = ap.parse_args()
    DATASET, MIN_BALANCE, MAX_STEPS, STATEFUL = args.dataset, args.min_balance, args.max_steps, args.stateful
    # post hoc conditions get their own output folders and run-id prefixes, so they never mix with the main runs
    suffix = ("_stateful" if STATEFUL else "") + (f"_{MAX_STEPS}steps" if MAX_STEPS != 15 else "")
    OUT = Path("results/task3" + ("_hao" if DATASET == "hao" else "") + suffix)
    RUNS_CSV = OUT / "task3_runs.csv"
    OUT.mkdir(parents=True, exist_ok=True)
    done = pd.read_csv(RUNS_CSV) if RUNS_CSV.exists() else pd.DataFrame(columns=["run_id", "model_key"])
    keys = args.models.split(",")
    # interleave models (run 1 of each, then run 2 of each, ...) so a budget stop leaves them balanced
    prefix = ("hao_" if DATASET == "hao" else "") + ("stateful_" if STATEFUL else "") + (f"s{MAX_STEPS}_" if MAX_STEPS != 15 else "")
    todo = [(k, f"{prefix}{k}_run{i:02d}") for i in range(1, args.runs + 1) for k in keys
            if f"{prefix}{k}_run{i:02d}" not in set(done.run_id)]  # full run id, prefix included
    print(f"{len(todo)} runs to do; balance ${balance():.2f}; prices $/M in/out:",
          {k: (round(v[0] * 1e6, 2), round(v[1] * 1e6, 2)) for k, v in PRICES.items()}, flush=True)
    for k, rid in todo:
        try:
            r = one_run(k, rid)
        except Stop as e:
            print(f"STOPPING ALL RUNS: {e}", flush=True)
            break
        print(f"{rid:<14} {r['status']:<6} steps {r['steps']:>2} errors {r['code_errors']} ${r['cost_usd']:.3f} "
              f"{r['seconds']:.0f}s | {r['stop_reason']}", flush=True)
        if r["stop_reason"].startswith("account balance would drop below"):
            print("STOPPING ALL RUNS: balance limit", flush=True)
            break
    print(f"balance now ${balance():.2f}")
