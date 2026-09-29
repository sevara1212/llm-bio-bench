"""Score Task 2: plain models, agent (if complete) and every baseline in one table.

Usage: python src/task2_score.py   (run task2_baselines.py first)

Metrics on the same 200 questions for every method:
  accuracy                    3-class (up / down / no_change)
  macro_F1                    over the 3 classes
  direction_acc_up_down       share of the 134 up/down questions predicted with the correct direction
  direction_acc_matched       same, on the 30 expression-matched up/down pairs (60 questions) from
                              data/task2_matched_pairs.json, where target expression level carries no
                              information about the direction
  said_no_change_on_up_down   share of up/down questions answered no_change
An answer that can't be parsed counts as wrong (it matches no class). Writes results/task2_scores.csv.
"""
import json
import os

import pandas as pd
from sklearn.metrics import f1_score

CLASSES = ["up", "down", "no_change"]
q = pd.DataFrame(json.load(open("data/task2_questions.json")))
pairs = json.load(open("data/task2_matched_pairs.json"))["pairs"]
matched = {p["up"] for p in pairs} | {p["down"] for p in pairs}

preds = {}
base = pd.read_csv("results/task2_baselines.csv").set_index("id")
for m in ["always_no_change", "coexpr_sign", "coexpr_sign_|r|>0.05", "collectri", "expression_only"]:
    name = "baseline: " + ("coexpr_sign_|r|>0.05 (fixed in advance)" if m == "coexpr_sign_|r|>0.05" else m)
    preds[name] = base[m]
MODELS = {"Claude Sonnet 5": "results/claude-sonnet-5_task2.csv",
          "GPT-5.6 Terra": "results/gpt-5.6-terra_task2.csv",
          "Gemini 3.8 Flash": "results/gemini-3.8-flash_task2.csv",
          "specialist agent (Claude + tools)": "results/agents_task2/agent-specialist_task2.csv",
          # post hoc condition A: the same agent with GPT / Gemini (OpenRouter, temperature 0)
          "specialist agent (GPT + tools, post hoc)": "results/agents_task2/agent-specialist-gpt-5.6-terra_task2.csv",
          "specialist agent (Gemini + tools, post hoc)": "results/agents_task2/agent-specialist-gemini-3.8-flash_task2.csv"}
unparsed = {}
for name, path in MODELS.items():
    if not os.path.exists(path):
        continue
    d = pd.read_csv(path)
    d = d[d.run == 0] if "run" in d else d
    if len(d) < len(q):
        print(f"skipping {name}: only {len(d)}/{len(q)} answers so far")
        continue
    d = d.set_index("id")
    preds[name] = d.predicted
    unparsed[name] = int(d.predicted.isna().sum())

rows = {}
qi = q.set_index("id")
ud = qi[qi.label != "no_change"]
mq = qi[qi.index.isin(matched)]
for name, p in preds.items():
    p = p.reindex(qi.index).fillna("unparsed")
    rows[name] = {
        "n": len(qi),
        "accuracy": (p == qi.label).mean(),
        "macro_F1": f1_score(qi.label, p, labels=CLASSES, average="macro", zero_division=0),
        "direction_acc_up_down": (p[ud.index] == ud.label).mean(),
        "direction_acc_matched": (p[mq.index] == mq.label).mean(),
        "said_no_change_on_up_down": (p[ud.index] == "no_change").mean(),
        "unparsed": unparsed.get(name, 0),
    }
t = pd.DataFrame(rows).T
pd.set_option("display.width", 200)
print(f"n = {len(qi)} questions ({(qi.label == 'up').sum()} up, {(qi.label == 'down').sum()} down, "
      f"{(qi.label == 'no_change').sum()} no_change); direction on {len(ud)} up/down, "
      f"matched pairs on {len(mq)} ({len(pairs)} pairs)\n")
print(t.round(3).to_string())
t.to_csv("results/task2_scores.csv")
