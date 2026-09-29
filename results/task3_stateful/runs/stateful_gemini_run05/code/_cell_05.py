# Check filtering thresholds
print("Cells with <200 genes:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Cells with >2500 genes:", (adata.obs['n_genes_by_counts'] > 2500).sum())
print("Cells with pct_counts_mt > 5%:", (adata.obs['pct_counts_mt'] > 5.0).sum())