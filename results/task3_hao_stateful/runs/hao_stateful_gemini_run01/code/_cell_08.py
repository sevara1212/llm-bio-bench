# Notice min n_genes is 499, max pct MT is 15.06%. All 7841 cells have already been pre-filtered for quality!
# Notice: min n_genes_by_counts is 499 (practically 500), max mt is 15%.
# This strongly suggests the dataset already had basic QC cell filtering applied!
# Let's verify genes:
sc.pp.filter_genes(adata, min_cells=3)
print("Genes remaining after min_cells=3:", adata.n_vars)