import scanpy as sc
adata_hvg = sc.read_h5ad('adata_hvg.h5ad')

sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden')

print(adata_hvg.obs['leiden'].value_counts())

adata_hvg.write('adata_clustered.h5ad')