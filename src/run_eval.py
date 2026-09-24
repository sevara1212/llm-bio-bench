"""Step 4: ask a model to name the cell type for each question.

Usage: python src/run_eval.py <model> [n_runs] [--dataset pbmc|hao]
  anthropic/claude-sonnet-5   -> Anthropic API directly (ANTHROPIC_API_KEY)
  openai/gpt-5.6-terra        -> OpenRouter (OPENROUTER_API_KEY)
  google/gemini-3.8-flash     -> OpenRouter

n_runs (default 1) asks every question that many times, to measure how much
answers wobble between identical calls. Logs tokens and cost per answer.
Resumable: (question, run) pairs already in the output CSV are skipped.
"""
import argparse
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

parser = argparse.ArgumentParser()
parser.add_argument("model", nargs="?", default="anthropic/claude-sonnet-5")
parser.add_argument("n_runs", nargs="?", type=int, default=1)
parser.add_argument("--dataset", default="pbmc", choices=["pbmc", "hao"])
args = parser.parse_args()
MODEL, N_RUNS, DATASET = args.model, args.n_runs, args.dataset
USE_CLAUDE_API = MODEL.startswith("anthropic/")
# OpenRouter caps new accounts at 20 requests/minute per model.
REQUESTS_PER_MIN = 50 if USE_CLAUDE_API else 18
WORKERS = 8
MAX_TOKENS = 4000  # same cap for every model; reasoning models need the room
# PBMC3k keeps its original wording so old and new runs stay comparable.
# The Hao prompt names the sample type, so answers like "neutrophil" count as ignoring context.
SAMPLE = {"pbmc": "human blood", "hao": "human peripheral blood mononuclear cells (PBMCs)"}[DATASET]
PROMPT = (
    f"These are the top marker genes of a cluster from {SAMPLE}. What cell type is it? "
    'Reply as JSON: {{"cell_type": ..., "confidence": 0-1}}\n\nGenes: {genes}'
)

if USE_CLAUDE_API:
    from claude_client import ask_claude
else:
    client = OpenAI(base_url="https://openrouter.ai/api/v1", api_key=os.environ["OPENROUTER_API_KEY"],
                    timeout=120, max_retries=6)


class RateLimiter:
    """Space requests out so at most `per_min` start in any minute."""

    def __init__(self, per_min):
        self.gap, self.next_at, self.lock = 60 / per_min, 0.0, threading.Lock()

    def wait(self):
        with self.lock:
            now = time.monotonic()
            start = max(now, self.next_at)
            self.next_at = start + self.gap
        time.sleep(start - now)


limiter = RateLimiter(REQUESTS_PER_MIN)


def call_openrouter(prompt):
    resp = client.chat.completions.create(
        model=MODEL,
        max_tokens=MAX_TOKENS,
        temperature=0,
        messages=[{"role": "user", "content": prompt}],
        extra_body={"usage": {"include": True}},  # OpenRouter adds cost to usage
    )
    usage = resp.usage.model_dump() if resp.usage else {}
    details = usage.get("completion_tokens_details") or {}
    return resp.choices[0].message.content or "", resp.choices[0].finish_reason, {
        "prompt_tokens": usage.get("prompt_tokens"),
        "completion_tokens": usage.get("completion_tokens"),
        "reasoning_tokens": details.get("reasoning_tokens"),
        "cost_usd": usage.get("cost"),
    }


def ask(q, run):
    limiter.wait()
    prompt = PROMPT.format(genes=", ".join(q["genes"]))
    if USE_CLAUDE_API:
        raw, finish_reason, usage = ask_claude(MODEL.split("/", 1)[1], prompt, MAX_TOKENS)
    else:
        raw, finish_reason, usage = call_openrouter(prompt)
    try:
        parsed = json.loads(re.search(r"\{.*\}", raw, re.S).group())
    except (AttributeError, json.JSONDecodeError):
        parsed = {}
    return {
        "id": q["id"], "format": q["format"], "knob": q["knob"], "level": q["level"],
        "replicate": q["replicate"], "run": run, "answer": q["answer"],
        "predicted": parsed.get("cell_type"), "confidence": parsed.get("confidence"),
        "raw": raw, "finish_reason": finish_reason, **usage,
    }


questions = json.load(open(f"data/questions_{DATASET}.json"))

# anthropic/claude-sonnet-5 --dataset hao -> results/claude-sonnet-5_hao.csv
out = f"results/{MODEL.split('/')[-1]}_{DATASET}.csv"
done = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame()
if "cost_usd" not in done.columns:  # old-format file from before cost logging: start over
    done = pd.DataFrame()
if not done.empty and "run" not in done.columns:  # files from before repeats were run 0
    done["run"] = 0
asked = set() if done.empty else set(zip(done["id"], done["run"]))
todo = [(q, run) for run in range(N_RUNS) for q in questions if (q["id"], run) not in asked]
print(f"{MODEL} on {DATASET} via {'Anthropic API' if USE_CLAUDE_API else 'OpenRouter'}: "
      f"{len(todo)} to ask, {len(done)} already done", flush=True)

CHECKPOINT_EVERY = 20  # save progress often, so a shutdown loses at most ~20 answers
rows = []


def save():
    results = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
    results.to_csv(out + ".tmp", index=False)
    os.replace(out + ".tmp", out)  # atomic: a crash mid-write can't corrupt the CSV
    return results


with ThreadPoolExecutor(WORKERS) as pool:
    futures = {pool.submit(ask, q, run): (q, run) for q, run in todo}
    for i, fut in enumerate(as_completed(futures), 1):
        q, run = futures[fut]
        try:
            rows.append(fut.result())
            print(f"[{i}/{len(todo)}] {q['id']:<40} pred={rows[-1]['predicted']}", flush=True)
            if len(rows) % CHECKPOINT_EVERY == 0:
                save()
        except Exception as e:  # leave it out; a rerun will retry it
            print(f"[{i}/{len(todo)}] {q['id']:<40} FAILED: {str(e)[:200]}", flush=True)

results = save()
print(f"\nSaved {len(results)}/{len(questions) * N_RUNS} answers to {out}  |  "
      f"total cost ${results.cost_usd.sum():.3f}  |  failed this run: {len(todo) - len(rows)}")
