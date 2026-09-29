# Let's inspect filtering
# Seurat standard tutorial uses:
# nFeature_RNA > 200 & nFeature_RNA < 2500
# percent.mt < 5
# Scanpy pbmc3k tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

print("Cells before filtering:", adata.n_obs)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
print("Cells after filtering:", adata.n_obs)
print("Genes after filtering:", adata.n_vars)