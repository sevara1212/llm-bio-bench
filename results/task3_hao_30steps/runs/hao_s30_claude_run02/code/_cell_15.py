import scanpy as sc
adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')
adata.obs['leiden'].value_counts().sort_index()