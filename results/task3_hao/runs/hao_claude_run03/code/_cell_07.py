import scanpy as sc

adata_hvg = sc.read_h5ad('pca.h5ad')

sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden')

print(adata_hvg.obs['leiden'].value_counts())

adata_hvg.write('clustered.h5ad')