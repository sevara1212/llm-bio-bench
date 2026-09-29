import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')

# Try lower resolution for cleaner clusters
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_05')
print(adata.obs['leiden_05'].value_counts())
adata.write('clustered.h5ad')