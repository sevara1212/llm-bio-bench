"""Agent conditions on the fixed Hao subset: Claude Sonnet 5 as a LangGraph prebuilt ReAct agent.

Usage: python src/run_agent.py <generic|specialist|specialist_nudge> [--limit N] [--questions subset|all_nonoise]

  generic     web search only (DuckDuckGo, via langchain-community)
  specialist  gene_info (MyGene.info: symbol or Ensembl ID -> symbol, name, summary)
              + the SAME CellMarker lookup as baseline_marker_lookup.py (Hao 2021 excluded, top 50
                markers per type, same mapping file): gene -> cell types, cell type -> markers
  specialist_nudge  identical to specialist plus a one-line system prompt (NUDGE below). Added after
                the error analysis of the specialist; the specialist itself is unchanged.

--questions subset (default) = the fixed 96-question subset; all_nonoise = all 600 non-noise Hao
questions (the set the plain models were compared on).

Everything else matches plain Claude: same model (claude-sonnet-5, Anthropic API), same prompt
(prompts.make_prompt), no temperature, default thinking, max_tokens 4000, same JSON answer parsing.
Limits: at most 5 tool calls per question (further calls get a "budget used up" reply) and a
2-minute wall-clock timeout per question.

Outputs:
  results/agents/<condition>/<question id>.json   full trajectory (every message, tool call, argument,
                                                   tool output, final answer, tokens, cost, latency)
  results/agents/agent-<condition>_hao.csv        one row per question, same columns as run_eval.py
Resumable: questions already in the CSV are skipped.
"""
import argparse
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout

import pandas as pd
import requests
from dotenv import load_dotenv
from langchain_anthropic import ChatAnthropic
from langchain_community.tools import DuckDuckGoSearchResults
from langchain_core.messages import AIMessage, HumanMessage, ToolMessage
from langchain_core.tools import StructuredTool
from langgraph.errors import GraphRecursionError
from langgraph.prebuilt import create_react_agent

from baseline_marker_lookup import load_markers
from prices import cost_usd
from prompts import make_prompt

load_dotenv()

parser = argparse.ArgumentParser()
parser.add_argument("condition", choices=["generic", "specialist", "specialist_nudge"])
parser.add_argument("--questions", choices=["subset", "all_nonoise"], default="subset")
parser.add_argument("--limit", type=int, help="pilot: only N questions, spread evenly over the subset")
args = parser.parse_args()

MODEL = "claude-sonnet-5"
DATASET = "hao"
MAX_TOOL_CALLS = 5
TIMEOUT_S = 120
WORKERS = 4
OUT_DIR = f"results/agents/{args.condition}"
CSV = f"results/agents/agent-{args.condition}_{DATASET}.csv"
os.makedirs(OUT_DIR, exist_ok=True)

# ---- shared data for the specialist tools -------------------------------------------------------
CM_MARKERS, _ = load_markers("cellmarker")  # {CellMarker cell name: top-50 symbols, strongest first}
CM_BY_LOWER = {ct.lower(): ct for ct in CM_MARKERS}
GENE_TO_TYPES = {}
for ct, ranked in CM_MARKERS.items():
    for rank, g in enumerate(ranked, start=1):
        GENE_TO_TYPES.setdefault(g, []).append((rank, ct))
# Same Ensembl -> symbol table the baseline uses.
_id_tables = [pd.read_csv(f"data/{DATASET}_markers.csv"), pd.read_csv(f"data/{DATASET}_background_genes.csv")]
ENSEMBL_TO_SYMBOL = {e: g for t in _id_tables for e, g in zip(t.ensembl, t.gene)}


def gene_info(genes: list[str]) -> str:
    """Look up genes on MyGene.info. Accepts gene symbols or Ensembl gene IDs (ENSG...), up to 20 at
    a time. Returns each gene's official symbol, full name and a short functional summary."""
    r = requests.post("https://mygene.info/v3/query", timeout=30, data={
        "q": ",".join(genes[:20]), "scopes": "symbol,ensembl.gene", "species": "human",
        "fields": "symbol,name,summary"})
    r.raise_for_status()
    out = {}
    for hit in r.json():  # MyGene can return several hits per query; prefer one with a summary
        q = hit["query"]
        if q in out and "notfound" not in out[q] and (out[q]["summary"] or not hit.get("summary")):
            continue
        out[q] = ({"notfound": True} if hit.get("notfound") else
                  {"symbol": hit.get("symbol"), "name": hit.get("name"),
                   "summary": (hit.get("summary") or "")[:400]})
    return json.dumps(out, separators=(",", ":"))


def cellmarker_gene_to_cell_types(genes: list[str]) -> str:
    """Look up genes in the CellMarker database (human blood, healthy donors). Accepts gene symbols or
    Ensembl gene IDs. For each gene, returns the cell types that list it among their top-50 markers,
    with the gene's rank in that cell type's list (1 = strongest marker). Up to 20 cell types per gene."""
    out = {}
    for g in genes[:20]:
        sym = str(ENSEMBL_TO_SYMBOL.get(g, g)).upper()
        hits = sorted(GENE_TO_TYPES.get(sym, []))
        out[g] = {"symbol": sym, "n_cell_types": len(hits),
                  "cell_types": [{"cell_type": ct, "rank": r} for r, ct in hits[:20]]}
    return json.dumps(out, separators=(",", ":"))


def cellmarker_cell_type_markers(cell_type: str) -> str:
    """Get the top-50 markers (strongest first) of one CellMarker cell type (human blood, healthy
    donors), e.g. "Naive CD4+ T cell". If the name isn't found exactly, returns similar names."""
    ct = CM_BY_LOWER.get(cell_type.strip().lower())
    if ct:
        return json.dumps({"cell_type": ct, "markers": CM_MARKERS[ct]})
    words = [w for w in re.split(r"\W+", cell_type.lower()) if len(w) > 2]
    similar = sorted(c for c in CM_MARKERS if any(w in c.lower() for w in words))[:25]
    return json.dumps({"not_found": cell_type, "similar_cell_types": similar})


ddg = DuckDuckGoSearchResults(output_format="list", num_results=5)


def web_search(query: str) -> str:
    """Search the web (DuckDuckGo). Returns the top 5 results with title, snippet and link."""
    return json.dumps(ddg.invoke(query), separators=(",", ":"))


SPECIALIST_TOOLS = [gene_info, cellmarker_gene_to_cell_types, cellmarker_cell_type_markers]
TOOL_FUNCS = {"generic": [web_search], "specialist": SPECIALIST_TOOLS, "specialist_nudge": SPECIALIST_TOOLS}
# The only difference between specialist and specialist_nudge. The other conditions have no system prompt.
NUDGE = "After converting any Ensembl IDs, always query CellMarker before answering."
SYSTEM_PROMPT = {"specialist_nudge": NUDGE}.get(args.condition)


def make_tools(counter):
    """Fresh tools per question, sharing one call counter so the 5-call budget covers all tools."""
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
            except Exception as e:  # report tool errors to the model instead of crashing the run
                return f"Tool error: {type(e).__name__}: {str(e)[:200]}"
        return StructuredTool.from_function(func=run, name=fn.__name__, description=fn.__doc__,
                                            args_schema=StructuredTool.from_function(fn).args_schema)
    return [wrap(f) for f in TOOL_FUNCS[args.condition]]


llm = ChatAnthropic(model=MODEL, max_tokens=4000, timeout=TIMEOUT_S, max_retries=6)


def serialize(m):
    d = {"type": m.type, "content": m.content}
    if isinstance(m, AIMessage):
        d["tool_calls"] = [{"name": c["name"], "args": c["args"], "id": c["id"]} for c in m.tool_calls]
        d["usage"] = m.usage_metadata
        d["stop_reason"] = (m.response_metadata or {}).get("stop_reason")
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
    agent = create_react_agent(llm, make_tools(counter), prompt=SYSTEM_PROMPT)
    prompt = make_prompt(DATASET, q["genes"])
    result, status = {"messages": [HumanMessage(prompt)]}, "ok"
    start = time.monotonic()
    pool = ThreadPoolExecutor(1)
    # recursion_limit is a safety net only; the tool budget above is the real limit.
    fut = pool.submit(agent.invoke, {"messages": [HumanMessage(prompt)]},
                      {"recursion_limit": 2 * MAX_TOOL_CALLS + 6})
    try:
        result = fut.result(timeout=TIMEOUT_S)
    except FuturesTimeout:
        status = "timeout"
    except GraphRecursionError:
        status = "recursion_limit"
    except Exception as e:  # e.g. network down after the client's retries: not saved, so a rerun retries it
        status = f"error: {type(e).__name__}: {str(e)[:120]}"
    pool.shutdown(wait=False)
    latency = time.monotonic() - start
    messages = result["messages"]
    raw = final_text(messages) if status == "ok" else ""
    try:
        parsed = json.loads(re.search(r"\{.*\}", raw, re.S).group())
    except (AttributeError, json.JSONDecodeError):
        parsed = {}
    ai = [m for m in messages if isinstance(m, AIMessage)]
    tin = sum((m.usage_metadata or {}).get("input_tokens", 0) for m in ai)
    tout = sum((m.usage_metadata or {}).get("output_tokens", 0) for m in ai)
    row = {
        "id": q["id"], "format": q["format"], "knob": q["knob"], "level": q["level"],
        "replicate": q["replicate"], "run": 0, "provider": "anthropic", "answer": q["answer"],
        "predicted": parsed.get("cell_type"), "confidence": parsed.get("confidence"), "raw": raw,
        "finish_reason": status, "prompt_tokens": tin, "completion_tokens": tout, "reasoning_tokens": None,
        "cost_usd": cost_usd(MODEL, tin, tout), "n_tool_calls": counter["n"],
        "n_model_calls": len(ai), "latency_s": round(latency, 1),
    }
    trajectory = {"condition": args.condition, "model": MODEL, "question": q, "system_prompt": SYSTEM_PROMPT,
                  "prompt": prompt,
                  "status": status, "messages": [serialize(m) for m in messages], "summary": row}
    json.dump(trajectory, open(f"{OUT_DIR}/{q['id']}.json", "w"), indent=1, default=str)
    return row


if args.questions == "subset":
    subset = json.load(open("data/agent_subset_hao.json"))["ids"]
else:
    subset = [q["id"] for q in json.load(open(f"data/questions_{DATASET}.json")) if q["knob"] != "noise"]
qs = {q["id"]: q for q in json.load(open(f"data/questions_{DATASET}.json"))}
todo_ids = subset[:: len(subset) // args.limit][: args.limit] if args.limit else subset
done = pd.read_csv(CSV) if os.path.exists(CSV) else pd.DataFrame()
todo = [qs[i] for i in todo_ids if done.empty or i not in set(done.id)]
print(f"agent-{args.condition}: {len(todo)} to ask, {len(done)} already done", flush=True)

rows, failed = [], 0
with ThreadPoolExecutor(WORKERS) as pool:
    for i, row in enumerate(pool.map(ask, todo), 1):
        if row["finish_reason"].startswith("error"):
            failed += 1
            print(f"[{i}/{len(todo)}] {row['id']:<45} FAILED (will retry on rerun): {row['finish_reason']}", flush=True)
            continue
        rows.append(row)
        print(f"[{i}/{len(todo)}] {row['id']:<45} tools={row['n_tool_calls']} {row['finish_reason']:<8} "
              f"{row['latency_s']:>5}s ${row['cost_usd']:.4f}  pred={row['predicted']}", flush=True)
        results = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
        results.to_csv(CSV + ".tmp", index=False)
        os.replace(CSV + ".tmp", CSV)
print(f"\nSaved {len(done) + len(rows)} rows to {CSV} | this run ${sum(r['cost_usd'] for r in rows):.3f} "
      f"| failed (not saved, rerun to retry): {failed}")
