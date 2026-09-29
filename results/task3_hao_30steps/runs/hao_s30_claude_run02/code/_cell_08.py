import scanpy as sc
adata = sc.read_h5ad('adata_filtered.h5ad')
adata_hvg = sc.read_h5ad('adata_hvg_clustered.h5ad')

# transfer clustering labels
adata.obs['leiden'] = adata_hvg.obs['leiden_05'].values
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = None

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.umap(adata_hvg)
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']

adata.write('adata_full_lognorm_clustered.h5ad')
print(adata)