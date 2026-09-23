"""Step 4: ask Claude to name the cell type for each question."""
import json
import re
import sys

import anthropic
import pandas as pd
from dotenv import load_dotenv

load_dotenv()

MODEL = sys.argv[1] if len(sys.argv) > 1 else "claude-sonnet-5"
PROMPT = (
    "These are the top marker genes of a cluster from human blood. What cell type is it? "
    'Reply as JSON: {{"cell_type": ..., "confidence": 0-1}}\n\nGenes: {genes}'
)

client = anthropic.Anthropic()
questions = json.load(open("data/questions.json"))

rows = []
for q in questions:
    msg = client.messages.create(
        model=MODEL,
        max_tokens=200,
        temperature=0,
        messages=[{"role": "user", "content": PROMPT.format(genes=", ".join(q["genes"]))}],
    )
    raw = msg.content[0].text
    try:
        parsed = json.loads(re.search(r"\{.*\}", raw, re.S).group())
    except (AttributeError, json.JSONDecodeError):
        parsed = {}
    rows.append({
        "id": q["id"], "format": q["format"], "answer": q["answer"],
        "predicted": parsed.get("cell_type"), "confidence": parsed.get("confidence"),
        "raw": raw,
    })
    print(f"{q['id']:<22} truth={q['answer']:<18} pred={parsed.get('cell_type')}")

out = f"results/{MODEL.split('-')[1] if MODEL.startswith('claude') else MODEL}_pbmc.csv"
out = "results/claude_pbmc.csv" if MODEL.startswith("claude") else f"results/{MODEL}_pbmc.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"\nSaved {len(rows)} answers to {out}")
