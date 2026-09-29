# Let's inspect min genes / counts
print("Min n_genes_by_counts:", adata.obs['n_genes_by_counts'].min())
print("Min total_counts:", adata.obs['total_counts'].min())
print("Genes expressed in < 3 cells:", (np.asarray((adata.X > 0).sum(axis=0)).squeeze() < 3).sum())