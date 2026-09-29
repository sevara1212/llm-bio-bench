# Let's inspect filtering criteria
# Typically in scanpy PBMC3k tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

print("Cells with n_genes < 200:", (adata.obs.n_genes_by_counts < 200).sum())
print("Cells with n_genes >= 2500:", (adata.obs.n_genes_by_counts >= 2500).sum())
print("Cells with pct_counts_mt >= 5:", (adata.obs.pct_counts_mt >= 5).sum())