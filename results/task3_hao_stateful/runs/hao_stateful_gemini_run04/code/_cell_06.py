# Check how many cells would be filtered by various criteria
print("Total cells:", adata.n_obs)
print("pct_counts_mt > 15:", (adata.obs['pct_counts_mt'] > 15).sum())
print("pct_counts_mt > 20:", (adata.obs['pct_counts_mt'] > 20).sum())
print("n_genes_by_counts < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("n_genes_by_counts < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())
print("n_genes_by_counts > 5000:", (adata.obs['n_genes_by_counts'] > 5000).sum())