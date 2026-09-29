import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print("Obs head:", adata.obs.head())
print("Var head:", adata.var.head())
print("X shape and type:", adata.X.shape, type(adata.X))
print("Any raw layers?", adata.layers.keys())