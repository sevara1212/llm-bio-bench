import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("X shape:", adata.shape)
print("obs head:\n", adata.obs.head())
print("var head:\n", adata.var.head())
print("X min/max/type:", adata.X.min(), adata.X.max(), type(adata.X))