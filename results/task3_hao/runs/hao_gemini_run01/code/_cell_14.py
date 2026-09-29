import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
# Save raw copy
raw_adata = adata.copy()

# Filter genes expressed in < 3 cells
sc.pp.filter_genes(adata, min_cells=3)
# Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Store normalized in raw for marker finding
adata.raw = adata

# HVG
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)

adata.obs['leiden'] = adata_hvg.obs['leiden']

# Check marker expression
markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'NCAM1', 'FCGR3A', 'CD14', 'MS4A1', 'CD79A', 'PPBP', 'MZB1', 'LILRA4', 'CLEC9A', 'CD1C', 'FCER1A', 'FOXP3', 'MKI67', 'TYMS', 'NKG7', 'GNLY', 'GZMB', 'GZMK', 'CCR7', 'IL7R', 'TCF7', 'SLC4A10', 'TRDC']

mean_exp = {}
for m in markers:
    if m in adata.var_names:
        mean_exp[m] = adata.obs_vector(m)
df_markers = pd.DataFrame(mean_exp, index=adata.obs_names)
df_markers['leiden'] = adata.obs['leiden']
grouped = df_markers.groupby('leiden').mean()
print(grouped[['CD3D', 'CD4', 'CD8A', 'MS4A1', 'CD14', 'FCGR3A', 'PPBP', 'NKG7', 'GNLY', 'MZB1', 'LILRA4', 'CLEC9A', 'CD1C', 'TYMS', 'MKI67']])