adata_proc = adata.copy()

# Store raw counts
adata_proc.raw = adata_proc

# Normalize
sc.pp.normalize_total(adata_proc, target_sum=1e4)
sc.pp.log1p(adata_proc)

# Highly variable genes
sc.pp.highly_variable_genes(adata_proc, n_top_genes=2000, flavor='seurat')

# PCA
sc.pp.pca(adata_proc, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state=42)

# Neighbors and clustering
sc.pp.neighbors(adata_proc, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata_proc, resolution=0.5, random_state=42, key_added='leiden_0.5')
sc.tl.leiden(adata_proc, resolution=0.8, random_state=42, key_added='leiden_0.8')

print("Cluster counts res 0.5:")
print(adata_proc.obs['leiden_0.5'].value_counts())
print("\nCluster counts res 0.8:")
print(adata_proc.obs['leiden_0.8'].value_counts())