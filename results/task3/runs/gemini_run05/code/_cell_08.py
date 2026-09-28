import scanpy as sc
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.filter_genes(adata, min_cells=3)
print("After filtering:", adata.shape)

# Normalize and log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store raw
adata.raw = adata

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes:", sum(adata.var.highly_variable))

# Scale and PCA
adata_processed = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_processed, max_value=10)
sc.tl.pca(adata_processed, svd_solver='arpack')

sc.pp.neighbors(adata_processed, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_processed)
sc.tl.leiden(adata_processed, resolution=0.5)

adata.obs['leiden_0.5'] = adata_processed.obs['leiden']
sc.tl.leiden(adata_processed, resolution=0.8)
adata.obs['leiden_0.8'] = adata_processed.obs['leiden']

print("Clusters at res 0.5:", adata.obs['leiden_0.5'].value_counts())
print("Clusters at res 0.8:", adata.obs['leiden_0.8'].value_counts())