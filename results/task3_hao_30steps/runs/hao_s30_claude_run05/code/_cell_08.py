import scanpy as sc

adata_hvg = sc.read_h5ad('clustered_hvg.h5ad')

# Try lower resolution for cleaner clusters
sc.tl.leiden(adata_hvg, resolution=0.5, key_added='leiden_05')
print(adata_hvg.obs['leiden_05'].value_counts())
adata_hvg.write('clustered_hvg.h5ad')