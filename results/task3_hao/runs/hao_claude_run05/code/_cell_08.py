import scanpy as sc

adata = sc.read_h5ad('adata_clustered.h5ad')

# Try lower resolution for canonical PBMC cell types
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_05')
print(adata.obs['leiden_05'].value_counts())
adata.write('adata_clustered.h5ad')