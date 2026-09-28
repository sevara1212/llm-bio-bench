import scanpy as sc
adata = sc.read_h5ad('processed.h5ad')

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=1.0, key_added='leiden', flavor='igraph', n_iterations=2)
sc.tl.umap(adata)

print(adata.obs['leiden'].value_counts())
adata.write('clustered.h5ad')