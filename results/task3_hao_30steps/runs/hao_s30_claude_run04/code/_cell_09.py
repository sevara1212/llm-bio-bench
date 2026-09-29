import scanpy as sc
import numpy as np
import pandas as pd

adata_hvg = sc.read_h5ad('adata_processed.h5ad')

# need raw log-normalized data with full genes, not just HVG subset, for marker genes
adata_full = sc.read_h5ad('adata_filtered.h5ad')
adata_full.layers['counts'] = adata_full.X.copy()
sc.pp.normalize_total(adata_full, target_sum=1e4)
sc.pp.log1p(adata_full)

adata_full.obs['leiden_r05'] = adata_hvg.obs['leiden_r05'].values
adata_full.obs['leiden_r03'] = adata_hvg.obs['leiden_r03'].values
adata_full.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata_full.obsm['X_umap'] = None

adata_full.write('adata_full_annotated.h5ad')
print(adata_full)