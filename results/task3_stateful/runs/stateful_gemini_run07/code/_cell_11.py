# Let's check marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_df = pd.DataFrame({group: [result['names'][group][i] for i in range(10)] for group in groups})
print("Top 10 markers per cluster:")
print(markers_df)