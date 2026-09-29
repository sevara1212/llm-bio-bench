print(adata.obs[adata.obs['pct_counts_hb'] > 5][['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_hb']].head(20))
print("Number of cells with pct_counts_hb > 5%:", (adata.obs['pct_counts_hb'] > 5).sum())
print("Number of cells with pct_counts_hb > 1%:", (adata.obs['pct_counts_hb'] > 1).sum())