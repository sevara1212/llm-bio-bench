import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
print(adata)
print(adata.X[:5,:5])
print(adata.obs.head())
print(adata.var.head())