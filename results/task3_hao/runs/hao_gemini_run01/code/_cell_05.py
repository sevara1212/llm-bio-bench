# Check distribution of genes, counts, mt pct
print("Genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Genes < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())
print("Genes > 5000:", (adata.obs['n_genes_by_counts'] > 5000).sum())
print("Counts > 40000:", (adata.obs['total_counts'] > 40000).sum())
print("pct_counts_mt > 10%:", (adata.obs['pct_counts_mt'] > 10).sum())
print("pct_counts_mt > 15%:", (adata.obs['pct_counts_mt'] > 15).sum())

# Gene detection frequency
n_cells_per_gene = np.sum(adata.X > 0, axis=0)
if hasattr(n_cells_per_gene, 'A1'):
    n_cells_per_gene = n_cells_per_gene.A1
print("Genes expressed in < 3 cells:", np.sum(n_cells_per_gene < 3))
print("Genes expressed in >= 3 cells:", np.sum(n_cells_per_gene >= 3))