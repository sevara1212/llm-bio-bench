import scanpy as sc
import pandas as pd
import numpy as np

# Load
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard scanpy pbmc3k tutorial filtering:
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

print("Before filtering:", adata.shape)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
print("After filtering:", adata.shape)

# Normalize and log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store raw/full log-normalized data for differential expression / plotting
adata.raw = adata

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes:", np.sum(adata.var.highly_variable))

# Regress out / scale
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)

# PCA
sc.tl.pca(adata_hvg, svd_solver='arpack')
sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_hvg)

# Let's try a few resolutions for leiden
for res in [0.4, 0.5, 0.6, 0.8]:
    sc.tl.leiden(adata_hvg, resolution=res, key_added=f'leiden_{res}')
    print(f"Res {res}: {len(adata_hvg.obs[f'leiden_{res}'].unique())} clusters")

adata_hvg.write('adata_processed.h5ad')