import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
print(adata.obs_keys(), adata.var_keys(), adata.uns_keys())