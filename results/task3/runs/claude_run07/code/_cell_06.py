import scanpy as sc
import numpy as np

adata = sc.read_h5ad('qc_filtered.h5ad')

# Save raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata  # store normalized log data for later marker gene reference

# HVG
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("n HVG:", adata.var.highly_variable.sum())

adata_hvg = adata[:, adata.var.highly_variable].copy()

# Scale
sc.pp.scale(adata_hvg, max_value=10)

# PCA
sc.tl.pca(adata_hvg, svd_solver='arpack')

adata_hvg.write('processed_hvg.h5ad')
adata.write('processed_full.h5ad')
print("done")