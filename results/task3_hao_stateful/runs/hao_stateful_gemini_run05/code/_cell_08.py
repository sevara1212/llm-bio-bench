# Find markers for clusters
# Let's test different resolutions, e.g. 0.3, 0.5
sc.tl.leiden(adata_hvg, resolution=0.3, key_added='leiden_0.3')
print("leiden 0.3 clusters:", adata_hvg.obs['leiden_0.3'].nunique())
print(adata_hvg.obs['leiden_0.3'].value_counts())

# Let's rank genes for leiden_0.5 and leiden_0.3
adata_proc.obs['leiden_0.3'] = adata_hvg.obs['leiden_0.3']
adata_proc.obs['leiden_0.5'] = adata_hvg.obs['leiden_0.5']

sc.tl.rank_genes_groups(adata_proc, groupby='leiden_0.3', method='wilcoxon')
markers_03 = pd.DataFrame(adata_proc.uns['rank_genes_groups']['names']).head(10)
print("Markers 0.3:\n", markers_03)