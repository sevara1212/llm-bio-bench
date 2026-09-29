import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("X min/max/type:", adata.X.min(), adata.X.max(), type(adata.X))
print("obs head:", adata.obs.head())
print("var head:", adata.var.head())