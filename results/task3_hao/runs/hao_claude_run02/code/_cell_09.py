import scanpy as sc
adata = sc.read_h5ad('clustered2.h5ad')
print(adata)
print(adata.obs['leiden'].value_counts())