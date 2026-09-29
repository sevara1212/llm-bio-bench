import scanpy as sc
import numpy as np

ad = sc.read_h5ad('normalized.h5ad')

ad_hvg = ad[:, ad.var.highly_variable].copy()

sc.pp.scale(ad_hvg, max_value=10)
sc.tl.pca(ad_hvg, n_comps=50, svd_solver='arpack')

# Keep PCA in the full object
ad.obsm['X_pca'] = ad_hvg.obsm['X_pca']

sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
sc.tl.leiden(ad, resolution=1.0, key_added='leiden')
sc.tl.umap(ad)

print(ad.obs['leiden'].value_counts())

ad.write('clustered.h5ad')