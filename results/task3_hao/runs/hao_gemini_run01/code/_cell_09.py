import scanpy as sc

adata = sc.read_h5ad('raw_counts.h5ad')
print("obs_names:", adata.obs_names[:5])
print("uns keys:", adata.uns.keys())
print("obsm keys:", adata.obsm.keys())