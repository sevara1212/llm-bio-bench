import scanpy as sc, os
adata=sc.read_h5ad('raw_counts.h5ad')
print(adata)
print('obs',adata.obs.head(), adata.obs.columns.tolist())
print('var',adata.var.head(), adata.var.columns.tolist())
print('X',adata.X.dtype, 'sum range',adata.X.sum(axis=1).min(),adata.X.sum(axis=1).max())
print('genes first',adata.var_names[:20].tolist())