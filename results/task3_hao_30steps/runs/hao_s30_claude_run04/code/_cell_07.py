import scanpy as sc
adata_hvg = sc.read_h5ad('adata_processed.h5ad')

# Try lower resolution
sc.tl.leiden(adata_hvg, resolution=0.5, key_added='leiden_r05')
print(adata_hvg.obs['leiden_r05'].value_counts())