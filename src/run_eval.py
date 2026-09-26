"""Step 4: ask a model to name the cell type for each question.

Usage: python src/run_eval.py <model> [n_runs] [--dataset pbmc|hao] [--skip-knob noise] [--test]
  anthropic/claude-sonnet-5   -> Anthropic API directly (ANTHROPIC_API_KEY)
  google/gemini-3.8-flash     -> Gemini API directly (GEMINI_API_KEY), OpenAI-compatible endpoint
  openai/gpt-5.6-terra        -> OpenRouter (OPENROUTER_API_KEY); so is anything else
  --provider openrouter       -> force OpenRouter, e.g. for Gemini while Google billing is down

--test asks the first question once, prints the raw response, tokens and cost, and saves nothing.

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

from prices import cost_usd
from prompts import make_prompt

load_dotenv()

parser = argparse.ArgumentParser()
parser.add_argument("model", nargs="?", default="anthropic/claude-sonnet-5")
parser.add_argument("n_runs", nargs="?", type=int, default=1)
parser.add_argument("--dataset", default="pbmc", choices=["pbmc", "hao"])
parser.add_argument("--skip-knob", action="append", default=[],
                    help="leave out a difficulty knob, e.g. --skip-knob noise (repeatable)")
parser.add_argument("--test", action="store_true", help="ask one question, print everything, save nothing")
parser.add_argument("--provider", choices=["anthropic", "google", "openrouter"],
                    help="override the provider picked from the model prefix, e.g. send google/... via openrouter")
args = parser.parse_args()
MODEL, N_RUNS, DATASET = args.model, args.n_runs, args.dataset
PROVIDER = args.provider or {"anthropic": "anthropic", "google": "google"}.get(MODEL.split("/")[0], "openrouter")
API_MODEL = MODEL if PROVIDER == "openrouter" else MODEL.split("/", 1)[1]  # "google/gemini-3.8-flash" -> "gemini-3.8-flash"
# Requests per minute, per provider. OpenRouter caps new accounts at 20/min per model.
# Gemini API limits depend on your billing tier - see aistudio.google.com/rate-limit - so 60 is a
# conservative start; the client also backs off and retries automatically on a 429.
REQUESTS_PER_MIN = {"anthropic": 50, "openrouter": 18, "google": 60}[PROVIDER]
WORKERS = 8
MAX_TOKENS = 4000  # same cap for every model; reasoning models need the room

if PROVIDER == "anthropic":
    from claude_client import ask_claude
elif PROVIDER == "google":
    client = OpenAI(base_url="https://generativelanguage.googleapis.com/v1beta/openai/",
                    api_key=os.environ["GEMINI_API_KEY"], timeout=120, max_retries=6)
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


def call_openai_compatible(prompt):
    """OpenRouter and the Gemini API both speak the OpenAI chat format."""
    resp = client.chat.completions.create(
        model=API_MODEL,
        max_tokens=MAX_TOKENS,
        temperature=0,
        messages=[{"role": "user", "content": prompt}],
        # OpenRouter adds the request's cost to usage; Google doesn't know this field.
        **({"extra_body": {"usage": {"include": True}}} if PROVIDER == "openrouter" else {}),
    )
    usage = resp.usage.model_dump() if resp.usage else {}
    details = usage.get("completion_tokens_details") or {}
    prompt_tokens = usage.get("prompt_tokens") or 0
    # Output incl. thinking = total - prompt. Google's completion_tokens may leave thinking out,
    # so derive it from the total, which always includes it (it's what Google bills as output).
    if usage.get("total_tokens"):
        output_tokens = usage["total_tokens"] - prompt_tokens
    else:
        output_tokens = usage.get("completion_tokens") or 0
    cost = usage.get("cost") if PROVIDER == "openrouter" else cost_usd(API_MODEL, prompt_tokens, output_tokens)
    return resp.choices[0].message.content or "", resp.choices[0].finish_reason, {
        "prompt_tokens": prompt_tokens,
        "completion_tokens": output_tokens,
        "reasoning_tokens": details.get("reasoning_tokens"),
        "cost_usd": cost,
    }, usage


def ask(q, run, test=False):
    limiter.wait()
    prompt = make_prompt(DATASET, q["genes"])
    if PROVIDER == "anthropic":
        raw, finish_reason, usage = ask_claude(API_MODEL, prompt, MAX_TOKENS)
        raw_usage = usage
    else:
        raw, finish_reason, usage, raw_usage = call_openai_compatible(prompt)
    if test:
        print(f"PROMPT:\n{prompt}\n\nRAW RESPONSE:\n{raw}\n\nFINISH REASON: {finish_reason}")
        print(f"\nUSAGE AS RETURNED BY {PROVIDER.upper()}:\n{json.dumps(raw_usage, indent=1, default=str)}")
        print(f"\nLOGGED: {json.dumps(usage)}")
    try:
        parsed = json.loads(re.search(r"\{.*\}", raw, re.S).group())
    except (AttributeError, json.JSONDecodeError):
        parsed = {}
    return {
        "id": q["id"], "format": q["format"], "knob": q["knob"], "level": q["level"],
        "replicate": q["replicate"], "run": run, "provider": PROVIDER, "answer": q["answer"],
        "predicted": parsed.get("cell_type"), "confidence": parsed.get("confidence"),
        "raw": raw, "finish_reason": finish_reason, **usage,
    }


questions = [q for q in json.load(open(f"data/questions_{DATASET}.json")) if q["knob"] not in args.skip_knob]

if args.test:
    print(f"TEST: {MODEL} -> provider {PROVIDER}, API model name {API_MODEL!r}, question {questions[0]['id']}\n")
    row = ask(questions[0], run=0, test=True)
    print(f"\nPARSED: cell_type={row['predicted']!r}, confidence={row['confidence']}  (truth: {row['answer']})")
    raise SystemExit

# anthropic/claude-sonnet-5 --dataset hao -> results/claude-sonnet-5_hao.csv
out = f"results/{MODEL.split('/')[-1]}_{DATASET}.csv"
done = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame()
if "cost_usd" not in done.columns:  # old-format file from before cost logging: start over
    done = pd.DataFrame()
if not done.empty and "run" not in done.columns:  # files from before repeats were run 0
    done["run"] = 0
if not done.empty and "provider" not in done.columns:  # rows from before this column existed
    done["provider"] = "anthropic" if PROVIDER == "anthropic" else "openrouter"
asked = set() if done.empty else set(zip(done["id"], done["run"]))
todo = [(q, run) for run in range(N_RUNS) for q in questions if (q["id"], run) not in asked]
print(f"{MODEL} on {DATASET} via {PROVIDER} (API model {API_MODEL!r}): "
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
