import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.raw = adata.copy()

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Normalization & log transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("HVGs:", adata.var['highly_variable'].sum())

# Keep raw normalized in raw or save
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added='leiden_0.5')
sc.tl.leiden(adata_hvg, resolution=0.8, key_added='leiden_0.8')
sc.tl.leiden(adata_hvg, resolution=1.0, key_added='leiden_1.0')

adata.obs['leiden_0.5'] = adata_hvg.obs['leiden_0.5']
adata.obs['leiden_0.8'] = adata_hvg.obs['leiden_0.8']
adata.obs['leiden_1.0'] = adata_hvg.obs['leiden_1.0']
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']

print("Leiden 0.5 clusters:", adata.obs['leiden_0.5'].value_counts())
print("Leiden 0.8 clusters:", adata.obs['leiden_0.8'].value_counts())
print("Leiden 1.0 clusters:", adata.obs['leiden_1.0'].value_counts())
adata.write_h5ad('processed.h5ad')