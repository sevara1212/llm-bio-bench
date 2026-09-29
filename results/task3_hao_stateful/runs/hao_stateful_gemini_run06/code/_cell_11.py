# Let's check marker genes for resolution 0.5
sc.tl.rank_genes_groups(adata_proc, groupby='leiden_0.5', method='wilcoxon')
markers = pd.DataFrame(adata_proc.uns['rank_genes_groups']['names']).head(10)
print(markers)