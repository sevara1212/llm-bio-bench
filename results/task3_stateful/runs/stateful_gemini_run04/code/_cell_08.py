# Let's find marker genes for leiden clusters
sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='t-test')
result = adata_filtered.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_df = pd.DataFrame({group: result['names'][group][:10] for group in groups})
print("Top 10 markers per cluster (res=0.5):")
print(markers_df)