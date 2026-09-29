import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("obs head:")
print(adata.obs.head())
print("var head:")
print(adata.var.head())
print("X dtype, shape, min, max, is sparse:", type(adata.X), adata.shape, adata.X[:10, :10].toarray() if hasattr(adata.X, 'toarray') else adata.X[:10, :10])