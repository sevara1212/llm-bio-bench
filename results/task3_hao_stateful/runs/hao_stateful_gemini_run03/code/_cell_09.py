# Let's check genes and cells filtering
# If we filter cells, how many cells would be filtered?
# Standard QC:
# min_genes: 200 (all cells have >= 499)
# pct_counts_mt < 20% (all cells have <= 15.06%)
# So if we use standard QC thresholds: min_genes=200 or 300, mt < 15 or 20, almost NO cells are removed!
# What about min_cells per gene? sc.pp.filter_genes(adata, min_cells=3)
# Let's check how many genes have < 3 cells.
genes_expressed = np.diff(adata.X.tocsc().indptr)
print("Genes with 0 cells:", (genes_expressed == 0).sum())
print("Genes with < 3 cells:", (genes_expressed < 3).sum())