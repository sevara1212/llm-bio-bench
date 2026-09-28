import scanpy as sc
import pandas as pd
import numpy as np

# Load
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
raw_barcodes = set(adata.obs_names)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Standard filtering as in Seurat / Scanpy PBMC3k tutorial
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.filter_genes(adata, min_cells=3)

# Normalize & log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# HVG
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# PCA
adata_proc = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_proc, max_value=10)
sc.tl.pca(adata_proc, svd_solver='arpack', random_state=0)

# Neighbors & clustering
sc.pp.neighbors(adata_proc, n_neighbors=10, n_pcs=40, random_state=0)
for res in [0.5, 0.6, 0.7, 0.8, 1.0]:
    sc.tl.leiden(adata_proc, resolution=res, key_added=f'leiden_{res}', random_state=0)
    adata.obs[f'leiden_{res}'] = adata_proc.obs[f'leiden_{res}']
    print(f"Res {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")

adata_proc.write('adata_proc.h5ad')
adata.write('adata_filtered.h5ad')