import scanpy as sc
import numpy as np

adata = sc.read_h5ad('normalized.h5ad')

adata_hvg = adata[:, adata.var['highly_variable']].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')

# copy PCA back to full adata
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')
print(adata.obs['leiden'].value_counts())

sc.tl.umap(adata)

adata.write('clustered.h5ad')