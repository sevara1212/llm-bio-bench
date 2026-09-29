print(adata.obs[['pct_counts_mt', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts']].describe())
print("\npct_counts_mt quantiles:")
print(adata.obs['pct_counts_mt'].quantile([0.01, 0.05, 0.5, 0.95, 0.99, 1.0]))