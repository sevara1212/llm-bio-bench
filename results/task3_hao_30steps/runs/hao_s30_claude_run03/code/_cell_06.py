import scanpy as sc
import numpy as np

ad = sc.read_h5ad('qc_filtered.h5ad')

# save raw counts layer
ad.layers['counts'] = ad.X.copy()

# Normalize
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
ad.raw = ad  # store normalized log data for later marker gene lookups

# HVG selection
sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor='seurat')
print(ad.var.highly_variable.sum())

ad.write('normalized.h5ad')