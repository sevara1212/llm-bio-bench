import scanpy as sc
import numpy as np

ad = sc.read_h5ad('qc_filtered.h5ad')

# save raw counts
ad.layers['counts'] = ad.X.copy()

# normalization
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
ad.raw = ad

# HVGs
sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor='seurat')
print(ad.var['highly_variable'].sum())

ad_hvg = ad[:, ad.var.highly_variable].copy()

sc.pp.scale(ad_hvg, max_value=10)
sc.tl.pca(ad_hvg, svd_solver='arpack', n_comps=50)

ad.obsm['X_pca'] = ad_hvg.obsm['X_pca']

sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
sc.tl.leiden(ad, resolution=1.0, key_added='leiden')

print(ad.obs['leiden'].value_counts())

ad.write('processed.h5ad')