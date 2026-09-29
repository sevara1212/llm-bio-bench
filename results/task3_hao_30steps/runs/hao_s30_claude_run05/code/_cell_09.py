import scanpy as sc

adata_full = sc.read_h5ad('processed_full.h5ad')
adata_hvg = sc.read_h5ad('clustered_hvg.h5ad')

adata_full.obs['leiden'] = adata_hvg.obs['leiden'].values
adata_full.obs['leiden_05'] = adata_hvg.obs['leiden_05'].values
adata_full.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata_full.obsm['X_umap'] = adata_hvg.obsm['X_umap']

adata_full.write('full_with_clusters.h5ad')
print(adata_full)