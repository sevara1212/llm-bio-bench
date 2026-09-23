"""Step 2: scanpy PBMC3k tutorial -> top markers per cluster + cluster labels."""
import pandas as pd
import scanpy as sc

N_TOP = 10

# Canonical markers from the scanpy PBMC3k tutorial, used to label clusters.
CELL_TYPE_MARKERS = {
    "CD4 T cells": ["IL7R", "CD3D", "LDHB"],
    "CD14 Monocytes": ["CD14", "LYZ", "S100A8"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "CD8 T cells": ["CD8A", "CD8B", "CCL5"],
    "NK cells": ["GNLY", "NKG7", "GZMB"],
    "FCGR3A Monocytes": ["FCGR3A", "MS4A7"],
    "Dendritic cells": ["FCER1A", "CST3"],
    "Megakaryocytes": ["PPBP", "PF4"],
}

adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
symbol_to_ensembl = adata.var["gene_ids"].to_dict()

# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5)].copy()

# Normalize, pick variable genes, scale
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# Cluster
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph", n_iterations=2, directed=False)

# Marker genes
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
names = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(N_TOP)
rows = []
for cluster in names.columns:
    for rank, gene in enumerate(names[cluster], start=1):
        rows.append({"cluster": cluster, "rank": rank, "gene": gene,
                     "ensembl": symbol_to_ensembl[gene]})
pd.DataFrame(rows).to_csv("data/pbmc_markers.csv", index=False)

# Label each cluster by which canonical marker set it expresses most.
# CHECK THIS BY EYE against the tutorial before trusting it.
expr = sc.get.obs_df(adata, keys=["leiden"] + sorted({g for gs in CELL_TYPE_MARKERS.values() for g in gs}),
                     use_raw=True).groupby("leiden", observed=True).mean()
expr = (expr - expr.mean()) / expr.std()  # compare each gene across clusters, not raw levels
labels = []
for cluster, row in expr.iterrows():
    scores = {ct: row[gs].mean() for ct, gs in CELL_TYPE_MARKERS.items()}
    labels.append({"cluster": cluster, "cell_type": max(scores, key=scores.get),
                   "n_cells": int((adata.obs.leiden == cluster).sum())})
labels = pd.DataFrame(labels)
labels.to_csv("data/pbmc_labels.csv", index=False)

print(labels.to_string(index=False))
print("\nTop 5 markers per cluster:")
print(names.head(5).to_string())
