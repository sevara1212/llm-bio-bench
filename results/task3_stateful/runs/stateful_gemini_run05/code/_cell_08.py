# Let's run the standard Scanpy PBMC 3k pipeline on adata_filtered:
adata_filtered.raw = adata_filtered

sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata_filtered.raw = adata_filtered # freeze normalized and log1p

adata_filtered_hvg = adata_filtered[:, adata_filtered.var.highly_variable].copy()
sc.pp.scale(adata_filtered_hvg, max_value=10)
sc.tl.pca(adata_filtered_hvg, svd_solver='arpack')
sc.pp.neighbors(adata_filtered_hvg, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_filtered_hvg)
sc.tl.leiden(adata_filtered_hvg)

adata_filtered.obs['leiden'] = adata_filtered_hvg.obs['leiden']
print(adata_filtered.obs['leiden'].value_counts())