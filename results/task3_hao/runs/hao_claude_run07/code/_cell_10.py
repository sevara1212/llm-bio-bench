import scanpy as sc
adata = sc.read_h5ad('adata_filtered.h5ad')
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')

adata_hvg2 = sc.read_h5ad('adata_hvg.h5ad')
adata.obsm['X_pca'] = adata_hvg2.obsm['X_pca']
adata.obs['leiden'] = adata_hvg2.obs['leiden_0.4'].values

print(adata.obs['leiden'].value_counts())
adata.write('adata_processed.h5ad')