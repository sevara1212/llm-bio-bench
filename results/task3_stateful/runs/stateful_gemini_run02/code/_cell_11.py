# Let's check cluster 0 vs 2, or subclustering of cluster 0 at higher resolution
# Notice cluster 0 is 1170 cells (T cells). Can we separate CD4 and CD8 T cells, or is resolution=1.0 or 1.2 better?
sc.tl.leiden(adata, resolution=1.0)
print(adata.obs['leiden'].value_counts())
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
top_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
print(top_genes)