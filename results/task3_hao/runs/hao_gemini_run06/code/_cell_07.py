import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Keep raw counts in adata.raw or layer
adata.layers['counts'] = adata.X.copy()

# Normalization & log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var['highly_variable'].sum())

# PCA
sc.tl.pca(adata, svd_solver='arpack')
# Neighbors & UMAP & Leiden
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden_0.5')
sc.tl.leiden(adata, resolution=0.8, random_state=42, key_added='leiden_0.8')

print("Leiden 0.5 clusters:", adata.obs['leiden_0.5'].value_counts())
print("Leiden 0.8 clusters:", adata.obs['leiden_0.8'].value_counts())