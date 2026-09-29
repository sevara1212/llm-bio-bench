import scanpy as sc, os
adata=sc.read_h5ad('raw_counts.h5ad')
print(adata)
print(adata.obs.head()); print(adata.var.head())
print(adata.X.min(),adata.X.max(),adata.X.sum())
print(adata.obs.columns.tolist(),adata.var.columns.tolist())