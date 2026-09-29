# Let's perform standard normalization, HVG, PCA, neighbors, clustering, and find markers!
# Let's write a complete pipeline to inspect the clusters and their markers.

adata.raw = adata.copy() # keep raw counts

# Normalization:
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
sc.pp.pca(adata, svd_solver='arpack', use_highly_variable=True, random_state=42)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)

print("Leiden clusters (res=0.5):")
print(adata.obs['leiden'].value_counts())

# Rank genes for groups
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
print("\nTop 10 marker genes per cluster:")
print(markers)