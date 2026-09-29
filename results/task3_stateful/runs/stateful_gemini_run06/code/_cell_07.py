# Check marker genes
sc.tl.rank_genes_groups(adata_scaled, 'leiden', method='t-test')
marker_genes = pd.DataFrame(adata_scaled.uns['rank_genes_groups']['names']).head(10)
print(marker_genes)