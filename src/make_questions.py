"""Step 3: questions per cluster x gene format x difficulty knob.

Usage: python src/make_questions.py [pbmc|hao]   (default pbmc)
  pbmc  PBMC3k, 8 Leiden clusters        -> data/questions_pbmc.json
  hao   Hao 2021 CITE-seq, 30 cell types -> data/questions_hao.json

Knobs (each question has `knob` and `level`):
  rank_window  which marker ranks to show: 1-10, 11-20, 21-30, 31-40
  n_genes      how many top genes to show: 1, 3, 5, 10
  noise        how many of the top 10 are swapped for random genes: 0, 2, 4, 6
  ribo_filter  top 10 with / without ribosomal + mitochondrial genes
"""
import json
import random
import re
import sys

import pandas as pd

NOISE_REPLICATES = 3  # different random draws per noise level, for error bars
RIBO_MT = re.compile(r"^(RPS|RPL|MT-)")

DATASET = sys.argv[1] if len(sys.argv) > 1 else "pbmc"

markers = pd.read_csv(f"data/{DATASET}_markers.csv", dtype={"cluster": str}).sort_values(["cluster", "rank"])
labels = pd.read_csv(f"data/{DATASET}_labels.csv", dtype={"cluster": str}).set_index("cluster")
background = pd.read_csv(f"data/{DATASET}_background_genes.csv")
bg_genes = list(zip(background.gene, background.ensembl))


def gene_sets(ranked):
    """Yield (knob, level, replicate, [(symbol, ensembl), ...]) for one cluster."""
    for start in (0, 10, 20, 30):
        yield "rank_window", f"{start + 1}-{start + 10}", 0, ranked[start:start + 10]
    for n in (1, 3, 5, 10):
        yield "n_genes", str(n), 0, ranked[:n]
    for n_noise in (0, 2, 4, 6):
        for rep in range(NOISE_REPLICATES if n_noise else 1):
            rng = random.Random(f"{ranked[0][0]}-{n_noise}-{rep}")
            genes = ranked[:10]
            for pos, fake in zip(rng.sample(range(10), n_noise), rng.sample(bg_genes, n_noise)):
                genes[pos] = fake
            yield "noise", str(n_noise), rep, genes
    yield "ribo_filter", "off", 0, ranked[:10]
    yield "ribo_filter", "on", 0, [g for g in ranked if not RIBO_MT.match(g[0])][:10]


questions = []
for cluster, group in markers.groupby("cluster"):
    ranked = list(zip(group.gene, group.ensembl))
    for knob, level, rep, genes in gene_sets(ranked):
        for fmt, idx in [("symbol", 0), ("ensembl", 1)]:
            questions.append({
                # "CD4 Naive" -> "CD4-Naive" so IDs have no spaces
                "id": f"{DATASET}_c{re.sub(r'[^A-Za-z0-9]+', '-', cluster)}_{knob}_{level}_r{rep}_{fmt}",
                "dataset": DATASET,
                "cluster": cluster,
                "genes": [g[idx] for g in genes],
                "format": fmt,
                "knob": knob,
                "level": level,
                "replicate": rep,
                "answer": labels.loc[cluster, "cell_type"],
            })

out = f"data/questions_{DATASET}.json"
with open(out, "w") as f:
    json.dump(questions, f, indent=2)
print(f"Wrote {len(questions)} questions to {out}")
print(pd.DataFrame(questions).groupby(["knob", "level"]).size().to_string())
