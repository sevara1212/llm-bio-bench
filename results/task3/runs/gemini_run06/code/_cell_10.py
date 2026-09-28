# Let's inspect the 8 clusters in Scanpy tutorial vs resolution 0.8 / 0.9 / 1.0
adata.obs['leiden_0.8'] = adata_sub.obs['leiden_0.8']
adata.obs['leiden_1.0'] = adata_sub.obs['leiden_1.0']

# Let's check marker genes expression across clusters for res 0.8 and 1.0
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'CST3', 'FCER1A', 'PPBP']

sc.tl.rank_genes_groups(adata, 'leiden_0.8', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_08 = pd.DataFrame({group: result['names'][group][:8] for group in groups})
print("Top 8 markers (res 0.8, wilcoxon):")
print(markers_08)

sc.tl.rank_genes_groups(adata, 'leiden_1.0', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_10 = pd.DataFrame({group: result['names'][group][:8] for group in groups})
print("\nTop 8 markers (res 1.0, wilcoxon):")
print(markers_10)