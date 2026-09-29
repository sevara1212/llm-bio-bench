import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("n_cells:", adata.n_obs)
print("pct_counts_mt > 15:", (adata.obs['pct_counts_mt'] > 15).sum())
print("pct_counts_mt > 10:", (adata.obs['pct_counts_mt'] > 10).sum())
print("pct_counts_mt > 20:", (adata.obs['pct_counts_mt'] > 20).sum())
print("n_genes_by_counts < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("n_genes_by_counts < 300:", (adata.obs['n_genes_by_counts'] < 300).sum())
print("n_genes_by_counts < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())
print("total_counts < 1000:", (adata.obs['total_counts'] < 1000).sum())