import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)

adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata_hvg, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)

adata.obs['leiden'] = adata_hvg.obs['leiden']
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.obsp['connectivities'] = adata_hvg.obsp['connectivities']
adata.obsp['distances'] = adata_hvg.obsp['distances']

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write_h5ad('processed.h5ad')
print("Processed and saved successfully!")