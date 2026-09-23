"""Step 4: ask a model (via OpenRouter) to name the cell type for each question.

Usage: python src/run_eval.py [openrouter-model-id]
  e.g. python src/run_eval.py anthropic/claude-sonnet-5
"""
import json
import os
import re
import sys

import pandas as pd
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

MODEL = sys.argv[1] if len(sys.argv) > 1 else "anthropic/claude-sonnet-5"
PROMPT = (
    "These are the top marker genes of a cluster from human blood. What cell type is it? "
    'Reply as JSON: {{"cell_type": ..., "confidence": 0-1}}\n\nGenes: {genes}'
)

client = OpenAI(base_url="https://openrouter.ai/api/v1", api_key=os.environ["OPENROUTER_API_KEY"],
                timeout=120, max_retries=2)
questions = json.load(open("data/questions.json"))

rows = []
for q in questions:
    resp = client.chat.completions.create(
        model=MODEL,
        max_tokens=4000,
        temperature=0,
        messages=[{"role": "user", "content": PROMPT.format(genes=", ".join(q["genes"]))}],
    )
    raw = resp.choices[0].message.content or ""
    try:
        parsed = json.loads(re.search(r"\{.*\}", raw, re.S).group())
    except (AttributeError, json.JSONDecodeError):
        parsed = {}
    rows.append({
        "id": q["id"], "format": q["format"], "answer": q["answer"],
        "predicted": parsed.get("cell_type"), "confidence": parsed.get("confidence"),
        "raw": raw, "finish_reason": resp.choices[0].finish_reason,
    })
    print(f"{q['id']:<22} truth={q['answer']:<18} pred={parsed.get('cell_type')}")

# anthropic/claude-sonnet-5 -> results/claude-sonnet-5_pbmc.csv
out = f"results/{MODEL.split('/')[-1]}_pbmc.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"\nSaved {len(rows)} answers to {out}")
