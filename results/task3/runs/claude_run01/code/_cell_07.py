import scanpy as sc
adata = sc.read_h5ad('processed.h5ad')

sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=1.0)

print(adata.obs['leiden'].value_counts())

adata.write('clustered.h5ad')