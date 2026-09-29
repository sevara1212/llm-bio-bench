import scanpy as sc
adata = sc.read_h5ad('normalized.h5ad')
print(adata.obs_names[:20].tolist())
# check prefixes
prefixes = set([x.split('_')[0] for x in adata.obs_names])
print(prefixes)