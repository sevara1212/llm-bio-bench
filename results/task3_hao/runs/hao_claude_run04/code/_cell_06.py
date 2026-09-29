import scanpy as sc
import numpy as np

adata = sc.read_h5ad('qc_filtered.h5ad')

# save raw counts layer
adata.layers['counts'] = adata.X.copy()

# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
print(adata.var['highly_variable'].sum())

adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')

sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden')

print(adata_hvg.obs['leiden'].value_counts())

# transfer to full adata
adata.obs['leiden'] = adata_hvg.obs['leiden']
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsp = adata_hvg.obsp
adata.uns['neighbors'] = adata_hvg.uns['neighbors']
adata.uns['pca'] = adata_hvg.uns['pca']

sc.tl.umap(adata)

adata.write('clustered.h5ad')