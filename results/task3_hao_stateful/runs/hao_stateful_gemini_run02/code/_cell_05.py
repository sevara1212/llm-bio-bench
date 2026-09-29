# Let's inspect where this dataset comes from or what barcodes they are
print("Sample of obs_names:", adata.obs_names[:20].tolist())
# Check gene filtering: genes with 0 counts
n_cells_per_gene = np.array((adata.X > 0).sum(axis=0)).ravel()
print(f"Genes with 0 counts: {(n_cells_per_gene == 0).sum()}")
print(f"Genes in < 3 cells: {(n_cells_per_gene < 3).sum()}")
print(f"Min n_genes_by_counts: {adata.obs['n_genes_by_counts'].min()}")