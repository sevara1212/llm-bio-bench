import scanpy as sc
import numpy as np

adata = sc.read_h5ad('qc_filtered.h5ad')
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # store lognorm full data for marker genes

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')

adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added='leiden', flavor='igraph', n_iterations=2, directed=False)

print(adata_hvg.obs['leiden'].value_counts())

adata.obs['leiden'] = adata_hvg.obs['leiden'].values
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
sc.tl.umap(adata_hvg)
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']

adata.write('clustered2.h5ad')