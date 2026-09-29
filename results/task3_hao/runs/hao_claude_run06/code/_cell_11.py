import scanpy as sc

adata = sc.read_h5ad('clustered.h5ad')
sc.tl.leiden(adata, resolution=0.2, key_added='leiden_final')
print(adata.obs['leiden_final'].value_counts())

adata.write('clustered_final.h5ad')