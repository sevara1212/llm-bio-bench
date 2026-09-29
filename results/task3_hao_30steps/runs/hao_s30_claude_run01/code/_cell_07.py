import scanpy as sc
import numpy as np

adata = sc.read_h5ad('qc_filtered.h5ad')

# Save raw counts layer
adata.layers['counts'] = adata.X.copy()

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# save normalized log data before HVG subsetting, for marker gene analysis later
adata.raw = adata

# HVG
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
print(adata.var.highly_variable.sum())

adata.write('normalized.h5ad')