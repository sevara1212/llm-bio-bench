import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Keep raw counts in a layer
adata.layers['counts'] = adata.X.copy()

# Normalization & log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVG count:", adata.var['highly_variable'].sum())

# Store normalized data in .raw
adata.raw = adata

# Scale
adata_pca = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_pca, max_value=10)

# PCA
sc.tl.pca(adata_pca, svd_solver='arpack', random_state=42)

# Neighbors and UMAP and Leiden
sc.pp.neighbors(adata_pca, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata_pca, random_state=42)
sc.tl.leiden(adata_pca, resolution=0.5, random_state=42)

adata.obs['leiden_0.5'] = adata_pca.obs['leiden']
adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
adata.obsm['X_umap'] = adata_pca.obsm['X_umap']
adata.uns['neighbors'] = adata_pca.uns['neighbors']

# Let's test different resolutions
sc.tl.leiden(adata, resolution=0.8, random_state=42, key_added='leiden_0.8')
sc.tl.leiden(adata, resolution=0.4, random_state=42, key_added='leiden_0.4')
sc.tl.leiden(adata, resolution=0.6, random_state=42, key_added='leiden_0.6')

print("Clusters at res 0.4:", adata.obs['leiden_0.4'].value_counts())
print("Clusters at res 0.5:", adata.obs['leiden_0.5'].value_counts())
print("Clusters at res 0.8:", adata.obs['leiden_0.8'].value_counts())