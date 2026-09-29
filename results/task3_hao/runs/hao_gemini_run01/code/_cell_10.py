import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering: min_genes=200, pct_counts_mt < 20 (or no cell filtering if all pass min_genes=200 and mt < 20%)
# Notice all cells have n_genes >= 499 and pct_counts_mt <= 15.06%.
# If we filter min_genes=200, min_cells=3:
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print("After filter:", adata.shape)

adata.raw = adata  # keep raw counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes:", np.sum(adata.var.highly_variable))

adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack', random_state=42)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata_hvg, random_state=42)
sc.tl.leiden(adata_hvg, resolution=0.5, random_state=42)

adata.obs['leiden_0.5'] = adata_hvg.obs['leiden']
print("Clusters at res 0.5:", adata.obs['leiden_0.5'].value_counts())

# Let's also check resolution 0.6 and 0.8
sc.tl.leiden(adata_hvg, resolution=0.8, random_state=42)
adata.obs['leiden_0.8'] = adata_hvg.obs['leiden']
print("Clusters at res 0.8:", adata.obs['leiden_0.8'].value_counts())