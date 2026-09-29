import scanpy as sc
adata_hvg = sc.read_h5ad('adata_processed.h5ad')

sc.tl.leiden(adata_hvg, resolution=0.3, key_added='leiden_r03')
print(adata_hvg.obs['leiden_r03'].value_counts())
adata_hvg.write('adata_processed.h5ad')