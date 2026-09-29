# Let's inspect Louvain clustering as in the original tutorial:
# In scanpy tutorial: sc.tl.louvain(adata)
sc.tl.louvain(adata_proc, random_state=0)
adata_filtered.obs['louvain'] = adata_proc.obs['louvain']
print("Louvain value counts:")
print(adata_filtered.obs['louvain'].value_counts())

sc.tl.rank_genes_groups(adata_filtered, 'louvain', method='t-test')
result_louvain = adata_filtered.uns['rank_genes_groups']
groups_l = result_louvain['names'].dtype.names
markers_df_l = pd.DataFrame({group: result_louvain['names'][group][:10] for group in groups_l})
print("\nTop 10 markers per louvain cluster:")
print(markers_df_l)