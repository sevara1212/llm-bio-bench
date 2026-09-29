import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)

adata.raw = adata

# Normalization & log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVG count:", adata.var['highly_variable'].sum())

# PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', random_state=42)

# Neighbors and clustering
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata, random_state=42)

# Test a couple resolutions
for res in [0.3, 0.4, 0.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}', random_state=42)
    print(f"Res {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")

# Let's inspect res=0.3 or 0.4