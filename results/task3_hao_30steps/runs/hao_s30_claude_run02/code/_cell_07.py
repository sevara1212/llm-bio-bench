import scanpy as sc
adata_hvg = sc.read_h5ad('adata_hvg_clustered.h5ad')

sc.tl.leiden(adata_hvg, resolution=0.5, key_added='leiden_05')
print(adata_hvg.obs['leiden_05'].value_counts())
adata_hvg.write('adata_hvg_clustered.h5ad')