"""Step 5: map free-text predictions to canonical cell types and score them.

Usage: python src/score.py [pbmc|hao] [--skip-knob noise]   -> results/scores_<dataset>.csv

pbmc: regex RULES below map answers to the 8 PBMC3k types.
hao:  the LLM judge (judge.py, run it first) maps answers to the 30 Hao types.

Each answer gets a score:
  1   = correct cell type
  0.5 = right lineage only - either too vague ("T cell" for CD4 T cells) or,
        on hao, the wrong subtype of the right lineage (CD4 TCM for CD4 Naive)
  0   = wrong
strict = share scoring 1; lenient = mean score.

Only (question, run) pairs answered by EVERY model are scored, so models are compared on the
same questions; the number kept is printed as n.
"""
import glob
import sys

import pandas as pd

args = [a for a in sys.argv[1:]]
SKIP_KNOBS = [args[i + 1] for i, a in enumerate(args) if a == "--skip-knob"]
positional = [a for i, a in enumerate(args) if a != "--skip-knob" and (i == 0 or args[i - 1] != "--skip-knob")]
DATASET = positional[0] if positional else "pbmc"

from scoring import scorer  # noqa: E402

normalize, score, vague_labels = scorer(DATASET)

frames = []
for path in sorted(glob.glob(f"results/*_{DATASET}.csv")):
    if path.split("/")[-1].startswith("scores_"):
        continue
    df = pd.read_csv(path)
    if "knob" not in df.columns:  # skip old-format files
        continue
    if "run" not in df.columns:
        df["run"] = 0
    df.insert(0, "model", path.split("/")[-1].removesuffix(f"_{DATASET}.csv"))
    frames.append(df)
scores = pd.concat(frames, ignore_index=True)
if SKIP_KNOBS:
    print(f"[{DATASET}] leaving out knob(s) {SKIP_KNOBS}: "
          f"{scores.knob.isin(SKIP_KNOBS).groupby(scores.model).sum().to_dict()} answers per model")
    scores = scores[~scores.knob.isin(SKIP_KNOBS)]

# Keep only (question, run) pairs every model answered.
n_models = scores.model.nunique()
per_pair = scores.groupby(["id", "run"]).model.nunique()
common = per_pair[per_pair == n_models].index
before = scores.groupby("model").size().to_dict()
scores = scores.set_index(["id", "run"]).loc[common].reset_index()
print(f"[{DATASET}] answers per model before filtering: {before}")
print(f"[{DATASET}] kept {len(common)} (question, run) pairs answered by all {n_models} models "
      f"-> n = {len(common)} per model\n")

scores["normalized"] = scores.predicted.map(normalize)
if (scores.normalized == "NOT JUDGED").any():
    sys.exit(f"{(scores.normalized == 'NOT JUDGED').sum()} answers not judged yet: run python src/judge.py")
scores["score"] = scores.apply(score, axis=1)
scores["strict"] = (scores.score == 1).astype(float)
scores["vague"] = scores.normalized.isin(vague_labels).astype(float)
if DATASET == "hao":
    from judge import rule_match
    scores["label_source"] = ["none" if not isinstance(p, str) else "rule" if rule_match(p) else "judge"
                              for p in scores.predicted]
scores.to_csv(f"results/scores_{DATASET}.csv", index=False)

pd.set_option("display.width", 160)
print(f"[{DATASET}] runs per question: {scores.groupby('model').run.nunique().to_dict()}\n")
print("Overall accuracy (strict = exact cell type; lenient = 0.5 credit for right lineage)\n")
table = scores.groupby(["model", "format"]).agg(n=("strict", "size"), strict=("strict", "mean"),
                                                lenient=("score", "mean")).round(3)
print(table.unstack().to_string())
if "label_source" in scores:
    print("\nHow answers were mapped to labels:\n")
    print(scores.groupby("model").label_source.value_counts().unstack(fill_value=0).to_string())
print("\nVague-answer rate (named a lineage but no subtype)\n")
print(scores.groupby("model").vague.mean().round(3).to_string())
print("\nCost and tokens per answer\n")
print(scores.groupby("model").agg(answers=("id", "size"), total_cost=("cost_usd", "sum"),
                                  cost_per_answer=("cost_usd", "mean"),
                                  output_tokens=("completion_tokens", "mean")).round(4).to_string())
if DATASET == "pbmc":
    print("\nUnmatched predictions (check whether RULES should cover them):")
    print(scores[scores.normalized == "other"].predicted.value_counts().head(15).to_string())
