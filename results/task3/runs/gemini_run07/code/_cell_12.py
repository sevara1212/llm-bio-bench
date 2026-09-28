import scanpy as sc
import pandas as pd
import numpy as np

# Load and process data
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Keep track of all raw barcodes
all_barcodes = adata.obs_names.tolist()

# Filtering
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.filter_genes(adata, min_cells=3)

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save raw
adata.raw = adata

# Filter HVG & Scale
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# PCA & Neighbors
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.6, random_state=0)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

result = adata.uns['rank_genes_groups']
for g in result['names'].dtype.names:
    print(f"Cluster {g} (n={(adata.obs['leiden']==g).sum()}):", list(result['names'][g][:8]))

# Check canonical markers expression per cluster
canonical_markers = ['IL7R', 'CCR7', 'S100A4', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'GNLY', 'NKG7', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
available_markers = [m for m in canonical_markers if m in adata.raw.var_names]
sc.tl.dendrogram(adata, groupby='leiden')
mean_expr = pd.DataFrame(index=adata.obs['leiden'].cat.categories)
for m in available_markers:
    mean_expr[m] = [adata.raw[:, m].X[adata.obs['leiden'] == c].mean() for c in adata.obs['leiden'].cat.categories]
print("\nMean expression of canonical markers:")
print(mean_expr.round(2))