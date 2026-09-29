# Normalize and find variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata  # store raw normalized counts

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes:", np.sum(adata.var.highly_variable))

adata_scaled = adata.copy()
sc.pp.scale(adata_scaled, max_value=10)
sc.tl.pca(adata_scaled, svd_solver='arpack')
sc.pp.neighbors(adata_scaled, n_neighbors=10, n_pcs=10)
sc.tl.leiden(adata_scaled, resolution=0.5)

print(adata_scaled.obs['leiden'].value_counts())