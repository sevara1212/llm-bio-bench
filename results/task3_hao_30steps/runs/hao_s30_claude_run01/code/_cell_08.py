import scanpy as sc
import numpy as np

adata = sc.read_h5ad('normalized.h5ad')

# Subset to HVGs for PCA/clustering
adata_hvg = adata[:, adata.var.highly_variable].copy()

# Scale
sc.pp.scale(adata_hvg, max_value=10)

# PCA
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')

# Neighbors
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)

# Leiden clustering
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden')

print(adata_hvg.obs.leiden.value_counts())

# UMAP for visualization sanity check
sc.tl.umap(adata_hvg)

adata_hvg.write('clustered_hvg.h5ad')