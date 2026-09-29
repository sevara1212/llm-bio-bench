print(adata.obs['pct_counts_mt'].describe())
print("pct_counts_mt quantiles:")
print(adata.obs['pct_counts_mt'].quantile([0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0]))