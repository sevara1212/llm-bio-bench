import scanpy as sc
adata=sc.read_h5ad('raw_counts.h5ad')
print(adata)
print(adata.obs.head()); print(adata.var.head())
print(adata.X[:3,:3])
print(adata.obs.columns.tolist(), adata.var.columns.tolist())
print(adata.n_obs,adata.n_vars)