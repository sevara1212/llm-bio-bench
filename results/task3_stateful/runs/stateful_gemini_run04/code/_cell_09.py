# Let's try louvain clustering as well or leiden at resolution=0.8
sc.tl.leiden(adata_proc, resolution=0.8, random_state=0)
adata_filtered.obs['leiden_08'] = adata_proc.obs['leiden']
print(adata_filtered.obs['leiden_08'].value_counts())
sc.tl.rank_genes_groups(adata_filtered, 'leiden_08', method='t-test')
result08 = adata_filtered.uns['rank_genes_groups']
groups08 = result08['names'].dtype.names
markers_df08 = pd.DataFrame({group: result08['names'][group][:10] for group in groups08})
print("\nTop 10 markers per cluster (res=0.8):")
print(markers_df08)