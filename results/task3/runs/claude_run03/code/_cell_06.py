import scanpy as sc
import numpy as np

adata = sc.read_h5ad('adata_hvg.h5ad')

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=1.0)

print(adata.obs['leiden'].value_counts())

adata.write('adata_clustered.h5ad')