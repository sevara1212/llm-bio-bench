import scanpy as sc
import numpy as np

adata = sc.read_h5ad('qc.h5ad')

# save raw counts
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # for later marker gene use with all genes

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var.highly_variable.sum())

adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)

sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)

adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=1.0)

print(adata.obs['leiden'].value_counts())

adata.write('clustered.h5ad')