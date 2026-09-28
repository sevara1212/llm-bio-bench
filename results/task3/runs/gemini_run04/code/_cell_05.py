import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()

# QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("n_genes_by_counts summary:")
print(adata.obs['n_genes_by_counts'].describe())
print("total_counts summary:")
print(adata.obs['total_counts'].describe())
print("pct_counts_mt summary:")
print(adata.obs['pct_counts_mt'].describe())

# Check how many cells would be filtered with standard PBMC3k thresholds:
# Seurat / standard scanpy tutorial:
# min_genes=200, max_genes=2500, pct_counts_mt < 5%
filt = (adata.obs['n_genes_by_counts'] >= 200) & (adata.obs['n_genes_by_counts'] < 2500) & (adata.obs['pct_counts_mt'] < 5)
print(f"Cells passing standard filter (200 <= n_genes < 2500, mt < 5%): {filt.sum()} / {len(adata)}")