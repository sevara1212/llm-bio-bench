import scanpy as sc
import numpy as np

adata = sc.read_h5ad('adata_filtered.h5ad')

# Save raw counts layer
adata.layers['counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
print(adata.var['highly_variable'].sum())

adata.raw = adata  # store full log-norm data before subsetting

adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')

adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

adata.write('adata_norm.h5ad')
print("done")