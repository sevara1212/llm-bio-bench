import scanpy as sc
import numpy as np

# Load full normalized data
adata_full = sc.read_h5ad('normalized.h5ad')
adata_clustered = sc.read_h5ad('clustered_hvg.h5ad')

# Transfer cluster labels
adata_full.obs['leiden'] = adata_clustered.obs['leiden'].values
adata_full.obsm['X_pca'] = adata_clustered.obsm['X_pca']
adata_full.obsm['X_umap'] = adata_clustered.obsm['X_umap']

adata_full.write('full_with_clusters.h5ad')
print(adata_full.obs.leiden.value_counts())