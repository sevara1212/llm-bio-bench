# Let's inspect standard scanpy PBMC tutorial:
# In scanpy PBMC3k tutorial:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
# How many cells remain?
filt = (adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5) & (adata.obs.n_genes_by_counts >= 200)
print("Remaining cells with standard tutorial filter:", filt.sum(), "out of", len(adata))
print("Excluded cells:", (~filt).sum())