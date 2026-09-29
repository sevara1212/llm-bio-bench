import scanpy as sc, os
adata=sc.read_h5ad('raw_counts.h5ad')
print(adata)
print('obs',adata.obs.head(),adata.obs.columns.tolist())
print('var',adata.var.head(),adata.var.columns.tolist())
print('X range',adata.X.min(),adata.X.max(), 'sum',adata.X.sum())
print('obs names',adata.obs_names[:5].tolist()); print('var names',adata.var_names[:10].tolist())