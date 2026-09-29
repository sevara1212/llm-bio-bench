# Let's inspect Seurat vs Scanpy tutorial filtering:
# In scanpy tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
# Let's see how many cells are filtered out:
print("Cells with n_genes >= 200:", (adata.obs.n_genes_by_counts >= 200).sum())
print("Cells with n_genes < 2500:", (adata.obs.n_genes_by_counts < 2500).sum())
print("Cells with mt < 5%:", (adata.obs.pct_counts_mt < 5).sum())
print("Cells with both (scanpy tutorial):", ((adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5)).sum())