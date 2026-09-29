import scanpy as sc
adata_hvg = sc.read_h5ad('adata_processed.h5ad')
print(adata_hvg.obs.columns.tolist())