"""Hao et al. 2021 CITE-seq PBMC reference -> top markers per fine cell type.

Unlike PBMC3k, we don't cluster: the authors' own labels (celltype.l2, 30 types
after dropping doublets) are the groups, so the answer key is expert-curated.
The file is opened read-only on disk ("backed") and only a subsample of up to
CELLS_PER_TYPE cells per type is loaded, so this fits in 8 GB of RAM.

Input:  data/hao2021_pbmc.h5ad  (CELLxGENE, 161,764 cells)
Output: data/hao_markers.csv, data/hao_labels.csv, data/hao_background_genes.csv
"""
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

CELLS_PER_TYPE = 500
N_TOP = 50
DROP = {"Doublet"}

full = ad.read_h5ad("data/hao2021_pbmc.h5ad", backed="r")
obs = full.obs[["celltype.l1", "celltype.l2", "celltype.l3"]].copy()
obs = obs[~obs["celltype.l2"].isin(DROP)]

rng = np.random.default_rng(0)
keep = []
for _, idx in obs.groupby("celltype.l2", observed=True).indices.items():
    keep.extend(rng.choice(idx, min(len(idx), CELLS_PER_TYPE), replace=False))
keep = np.sort(full.obs_names.get_indexer(obs.index[keep]))  # positions within `obs` -> within the file

# Raw counts for the subsample only (CELLxGENE keeps counts in .raw, normalized values in .X).
counts = full.raw.X[keep] if full.raw is not None else full.X[keep]
var = full.var[["feature_name", "feature_is_filtered"]]
adata = ad.AnnData(X=counts, obs=full.obs.iloc[keep][["celltype.l1", "celltype.l2", "celltype.l3"]].copy(),
                   var=pd.DataFrame({"ensembl": var.index}, index=var.feature_name.astype(str).values))
adata = adata[:, ~var.feature_is_filtered.values].copy()
adata.var_names_make_unique()
full.file.close()
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes, {adata.obs['celltype.l2'].nunique()} cell types")

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, "celltype.l2", method="wilcoxon")

names = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(N_TOP)
symbol_to_ensembl = adata.var["ensembl"].to_dict()
rows = [{"cluster": ct, "rank": r, "gene": g, "ensembl": symbol_to_ensembl[g]}
        for ct in names.columns for r, g in enumerate(names[ct], start=1)]
pd.DataFrame(rows).to_csv("data/hao_markers.csv", index=False)

labels = (adata.obs.groupby("celltype.l2", observed=True)
          .agg(l1=("celltype.l1", "first"), n_cells=("celltype.l1", "size")).reset_index()
          .rename(columns={"celltype.l2": "cluster"}))
labels["cell_type"] = labels["cluster"]
labels["n_cells_total"] = labels.cluster.map(obs["celltype.l2"].value_counts())
labels[["cluster", "cell_type", "l1", "n_cells", "n_cells_total"]].to_csv("data/hao_labels.csv", index=False)

top_genes = {r["gene"] for r in rows}
background = [g for g in adata.var_names if g not in top_genes]
pd.DataFrame({"gene": background, "ensembl": [symbol_to_ensembl[g] for g in background]}
             ).to_csv("data/hao_background_genes.csv", index=False)

print(labels.to_string(index=False))
print("\nTop 5 markers per type:")
print(names.head(5).T.to_string())
