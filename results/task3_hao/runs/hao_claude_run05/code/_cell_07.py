import scanpy as sc

adata = sc.read_h5ad('adata_norm.h5ad')

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, use_rep='X_pca')
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')
sc.tl.umap(adata)

print(adata.obs['leiden'].value_counts())

adata.write('adata_clustered.h5ad')