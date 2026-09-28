import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("Summary of n_genes_by_counts:")
print(adata.obs['n_genes_by_counts'].describe())
print("\nSummary of total_counts:")
print(adata.obs['total_counts'].describe())
print("\nSummary of pct_counts_mt:")
print(adata.obs['pct_counts_mt'].describe())