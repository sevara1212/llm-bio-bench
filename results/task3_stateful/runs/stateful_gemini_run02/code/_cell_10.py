# In standard scanpy tutorial, resolution=0.8 or 0.6 is used to split CD4 T, CD8 T, monocytes, etc.
# Let's test resolution=0.8
sc.tl.leiden(adata, resolution=0.8)
print("Clusters at res 0.8:")
print(adata.obs['leiden'].value_counts())

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
top_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
print(top_genes)