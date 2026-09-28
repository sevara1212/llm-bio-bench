import scanpy as sc

adata_hvg = sc.read_h5ad('processed_hvg.h5ad')

sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata_hvg, resolution=1.0, flavor='igraph', n_iterations=2)

print(adata_hvg.obs.leiden.value_counts())

adata_hvg.write('clustered.h5ad')