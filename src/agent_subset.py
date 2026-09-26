"""Fixed random subset of Hao questions for the agent conditions.

Design (seed 0, so it is reproducible):
  - 48 base questions, each asked in BOTH formats (symbol + Ensembl) -> 96 questions, paired by format.
  - 12 base questions per difficulty knob, spread evenly over its levels:
      rank_window 3 x 4 levels, n_genes 3 x 4 levels, noise 3 x 4 levels, ribo_filter 6 x 2 levels.
  - Cell types are dealt out from a shuffled list of all 30, so each appears once or twice.
Noise is included: plain Claude and the lookup baseline answered all 1,200 Hao questions.

Writes data/agent_subset_hao.json (question IDs + design) - every agent condition uses exactly these.
"""
import json
import random
from collections import defaultdict

SEED = 0
PER_LEVEL = {"rank_window": 3, "n_genes": 3, "noise": 3, "ribo_filter": 6}

rng = random.Random(SEED)
questions = json.load(open("data/questions_hao.json"))
by_key = defaultdict(list)  # (knob, level, cluster) -> symbol-format questions (one per noise replicate)
for q in questions:
    if q["format"] == "symbol":
        by_key[(q["knob"], q["level"], q["cluster"])].append(q)
ensembl_of = {q["id"].removesuffix("_ensembl"): q["id"] for q in questions if q["format"] == "ensembl"}

clusters = sorted({q["cluster"] for q in questions})
deck = []
chosen = []
for knob in PER_LEVEL:
    levels = sorted({q["level"] for q in questions if q["knob"] == knob}, key=lambda l: (len(l), l))
    for level in levels:
        for _ in range(PER_LEVEL[knob]):
            if not deck:
                deck = clusters[:]
                rng.shuffle(deck)
            cluster = deck.pop()
            q = rng.choice(by_key[(knob, level, cluster)])  # picks a random noise replicate where there are 3
            chosen += [q["id"], ensembl_of[q["id"].removesuffix("_symbol")]]

json.dump({"seed": SEED, "design": PER_LEVEL, "n": len(chosen), "ids": chosen},
          open("data/agent_subset_hao.json", "w"), indent=1)
qs = {q["id"]: q for q in questions}
sub = [qs[i] for i in chosen]
print(f"{len(chosen)} questions saved to data/agent_subset_hao.json")
counts = defaultdict(int)
for q in sub:
    counts[(q["knob"], q["format"])] += 1
print("per knob x format:", dict(counts))
types = defaultdict(int)
for q in sub:
    types[q["answer"]] += 1
print(f"cell types covered: {len(types)}/30, questions per type: min {min(types.values())}, max {max(types.values())}")
