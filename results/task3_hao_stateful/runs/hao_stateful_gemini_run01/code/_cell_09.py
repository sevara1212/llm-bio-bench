adata.raw = adata.copy()

# Normalization and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var['highly_variable'].sum())

# Store normalized counts in raw
adata.raw = adata

# PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', random_state=42)

# Neighbors and clustering
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden_0.5')
sc.tl.leiden(adata, resolution=0.8, random_state=42, key_added='leiden_0.8')

print("Leiden 0.5 clusters:", adata.obs['leiden_0.5'].value_counts())
print("Leiden 0.8 clusters:", adata.obs['leiden_0.8'].value_counts())