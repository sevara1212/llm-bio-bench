import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()

# Check QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("n_genes_by_counts summary:")
print(adata.obs['n_genes_by_counts'].describe())
print("\ntotal_counts summary:")
print(adata.obs['total_counts'].describe())
print("\npct_counts_mt summary:")
print(adata.obs['pct_counts_mt'].describe())

print("\nCells with < 200 genes:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Cells with > 2500 genes:", (adata.obs['n_genes_by_counts'] > 2500).sum())
print("Cells with > 5% MT:", (adata.obs['pct_counts_mt'] > 5).sum())