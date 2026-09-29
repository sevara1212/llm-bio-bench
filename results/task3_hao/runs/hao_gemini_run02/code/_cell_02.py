import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("X type:", type(adata.X))
print("obs head:\n", adata.obs.head())
print("var head:\n", adata.var.head())
print("obs_names:", adata.obs_names[:5])
print("var_names:", adata.var_names[:5])
if adata.X.shape[0] > 0:
    print("X min/max:", adata.X.min(), adata.X.max())