"""Task 3 on Hao 2021: raw-count subset for the agent + the answer key (kept outside the agent's folder).

Usage: python src/task3_hao_data.py
  data/task3_hao_raw_counts.h5ad       what the agent gets: raw RNA counts only (X, integers), gene symbols as
                                       var_names plus gene_ids; cell barcodes; NO labels, NO protein (ADT) data,
                                       NO embeddings or other metadata. Cells shuffled (seed 0) so their order
                                       carries no information about type. (gitignored)
  data/task3_hao_expert_cells.csv      answer key: barcode, celltype.l1, celltype.l2 (never copied to the agent)
Up to 300 cells per celltype.l2 type (30 types; doublets dropped), seed 0.
"""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

CELLS_PER_TYPE, SEED = 300, 0
full = ad.read_h5ad("data/hao2021_pbmc.h5ad", backed="r")
obs = full.obs[["celltype.l1", "celltype.l2"]]
obs = obs[obs["celltype.l2"] != "Doublet"]
rng = np.random.default_rng(SEED)
keep = []
for _, idx in obs.groupby("celltype.l2", observed=True).indices.items():
    keep.extend(rng.choice(idx, min(len(idx), CELLS_PER_TYPE), replace=False))
keep = full.obs_names.get_indexer(obs.index[np.array(keep)])
keep = rng.permutation(keep)                       # shuffled order
order = np.argsort(keep)                           # backed reads need sorted indices
counts = sp.csr_matrix(full.raw.X[np.sort(keep)])[np.argsort(order)].astype(np.int32)
barcodes = full.obs_names[keep]
symbols = pd.Index(full.var.feature_name.astype(str).values)
var = pd.DataFrame({"gene_ids": full.var.index.values}, index=symbols)
var.index = ad.utils.make_index_unique(var.index)
agent = ad.AnnData(X=counts, obs=pd.DataFrame(index=barcodes), var=var)
agent.write_h5ad("data/task3_hao_raw_counts.h5ad", compression="gzip")
key = pd.DataFrame({"barcode": barcodes, "l1": full.obs["celltype.l1"].values[keep].astype(str),
                    "l2": full.obs["celltype.l2"].values[keep].astype(str)})
key.to_csv("data/task3_hao_expert_cells.csv", index=False)
print(f"{agent.n_obs:,} cells x {agent.n_vars:,} genes; {key.l2.nunique()} types; cells per type "
      f"{key.l2.value_counts().min()}-{key.l2.value_counts().max()}")
print("agent file contents -> obs columns:", list(agent.obs.columns), "| obsm:", list(agent.obsm), "| layers:",
      list(agent.layers), "| uns:", list(agent.uns), "| var columns:", list(agent.var.columns))
