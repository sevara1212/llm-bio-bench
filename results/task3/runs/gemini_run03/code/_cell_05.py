# Standard filtering: min_genes=200, min_cells=3, pct_counts_mt < 5% (standard Seurat/scanpy tutorial)
# Let's check how many cells pass:
print("Cells with n_genes >= 200:", (adata.obs['n_genes_by_counts'] >= 200).sum())
print("Cells with n_genes < 2500:", (adata.obs['n_genes_by_counts'] < 2500).sum())
print("Cells with pct_counts_mt < 5%:", (adata.obs['pct_counts_mt'] < 5).sum())
print("Passing all three:", ((adata.obs['n_genes_by_counts'] >= 200) & (adata.obs['n_genes_by_counts'] < 2500) & (adata.obs['pct_counts_mt'] < 5)).sum())