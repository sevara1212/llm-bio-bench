"""Step 3: one question per cluster per gene format (symbol / ensembl)."""
import json

import pandas as pd

markers = pd.read_csv("data/pbmc_markers.csv", dtype={"cluster": str})
labels = pd.read_csv("data/pbmc_labels.csv", dtype={"cluster": str}).set_index("cluster")

questions = []
for cluster, group in markers.sort_values("rank").groupby("cluster", sort=False):
    for fmt, col in [("symbol", "gene"), ("ensembl", "ensembl")]:
        questions.append({
            "id": f"pbmc_c{cluster}_{fmt}",
            "cluster": cluster,
            "genes": group[col].tolist(),
            "format": fmt,
            "answer": labels.loc[cluster, "cell_type"],
        })

with open("data/questions.json", "w") as f:
    json.dump(questions, f, indent=2)
print(f"Wrote {len(questions)} questions to data/questions.json")
