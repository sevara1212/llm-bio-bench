"""Step 4: ask a model to name the cell type for each question.

Usage: python src/run_eval.py <model>
  anthropic/claude-sonnet-5   -> Anthropic API directly (ANTHROPIC_API_KEY)
  openai/gpt-5.6-terra        -> OpenRouter (OPENROUTER_API_KEY)
  google/gemini-3.8-flash     -> OpenRouter

Logs tokens and cost per answer. Resumable: questions already in the output
CSV are skipped, so rerunning after a crash only asks the missing ones.
"""
import json
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

MODEL = sys.argv[1] if len(sys.argv) > 1 else "anthropic/claude-sonnet-5"
USE_CLAUDE_API = MODEL.startswith("anthropic/")
# OpenRouter caps new accounts at 20 requests/minute per model.
REQUESTS_PER_MIN = 50 if USE_CLAUDE_API else 18
WORKERS = 8
MAX_TOKENS = 4000  # same cap for every model; reasoning models need the room
PROMPT = (
    "These are the top marker genes of a cluster from human blood. What cell type is it? "
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


def ask(q):
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
        "replicate": q["replicate"], "answer": q["answer"],
        "predicted": parsed.get("cell_type"), "confidence": parsed.get("confidence"),
        "raw": raw, "finish_reason": finish_reason, **usage,
    }


questions = json.load(open("data/questions.json"))

# anthropic/claude-sonnet-5 -> results/claude-sonnet-5_pbmc.csv
out = f"results/{MODEL.split('/')[-1]}_pbmc.csv"
done = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame()
if "cost_usd" not in done.columns:  # old-format file from before cost logging: start over
    done = pd.DataFrame()
todo = [q for q in questions if done.empty or q["id"] not in set(done["id"])]
print(f"{MODEL} via {'Anthropic API' if USE_CLAUDE_API else 'OpenRouter'}: "
      f"{len(todo)} to ask, {len(done)} already done", flush=True)

rows = []
with ThreadPoolExecutor(WORKERS) as pool:
    futures = {pool.submit(ask, q): q for q in todo}
    for i, fut in enumerate(as_completed(futures), 1):
        q = futures[fut]
        try:
            rows.append(fut.result())
            print(f"[{i}/{len(todo)}] {q['id']:<40} pred={rows[-1]['predicted']}", flush=True)
        except Exception as e:  # leave it out; a rerun will retry it
            print(f"[{i}/{len(todo)}] {q['id']:<40} FAILED: {str(e)[:200]}", flush=True)

results = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
results.to_csv(out, index=False)
print(f"\nSaved {len(results)}/{len(questions)} answers to {out}  |  "
      f"total cost ${results.cost_usd.sum():.3f}  |  failed this run: {len(todo) - len(rows)}")
