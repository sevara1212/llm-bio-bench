import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
adata.var_names_make_unique()

# Quality control
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(adata.obs.describe())
print(f"Cells with mt > 5%: {(adata.obs['pct_counts_mt'] > 5).sum()}")
print(f"Cells with n_genes_by_counts < 200: {(adata.obs['n_genes_by_counts'] < 200).sum()}")
print(f"Cells with n_genes_by_counts > 2500: {(adata.obs['n_genes_by_counts'] > 2500).sum()}")