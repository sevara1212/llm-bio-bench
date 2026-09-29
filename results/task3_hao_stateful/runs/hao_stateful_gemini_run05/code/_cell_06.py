# Let's save the raw counts in raw or keep adata intact
adata_proc = adata.copy()
sc.pp.filter_genes(adata_proc, min_cells=3)
sc.pp.normalize_total(adata_proc, target_sum=1e4)
sc.pp.log1p(adata_proc)
adata_proc.raw = adata_proc

sc.pp.highly_variable_genes(adata_proc, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", np.sum(adata_proc.var.highly_variable))

adata_hvg = adata_proc[:, adata_proc.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added='leiden_0.5')
sc.tl.leiden(adata_hvg, resolution=0.8, key_added='leiden_0.8')
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden_1.0')

print(adata_hvg.obs['leiden_0.5'].value_counts())
print(adata_hvg.obs['leiden_0.8'].value_counts())