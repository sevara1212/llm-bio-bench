import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("X type:", type(adata.X))
print("obs head:", adata.obs.head())
print("var head:", adata.var.head())
print("obs index name:", adata.obs_names[:5])
print("var index name:", adata.var_names[:5])