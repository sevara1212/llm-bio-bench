import scanpy as sc
adata = sc.read_h5ad('adata_filtered.h5ad')

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
print(adata.var['highly_variable'].sum())

adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)

sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden')
print(adata_hvg.obs['leiden'].value_counts())

adata_hvg.write('adata_hvg_clustered.h5ad')