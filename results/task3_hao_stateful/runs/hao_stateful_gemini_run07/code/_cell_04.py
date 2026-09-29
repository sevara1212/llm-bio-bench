# Look at distributions of QC metrics
print("pct_counts_mt quantiles:")
print(adata.obs['pct_counts_mt'].quantile([0.5, 0.9, 0.95, 0.98, 0.99, 1.0]))

print("\nn_genes_by_counts quantiles:")
print(adata.obs['n_genes_by_counts'].quantile([0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999]))

print("\ntotal_counts quantiles:")
print(adata.obs['total_counts'].quantile([0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999]))

print("\npct_counts_hb quantiles:")
print(adata.obs['pct_counts_hb'].quantile([0.9, 0.95, 0.98, 0.99, 1.0]))