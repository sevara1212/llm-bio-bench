import scanpy as sc
import numpy as np
adata = sc.read_h5ad('raw_counts.h5ad')

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], inplace=True, percent_top=None)

# basic filters
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata.shape)

# mt filter
print((adata.obs['pct_counts_mt']>15).sum())
print((adata.obs['pct_counts_mt']>10).sum())
print(adata.obs['total_counts'].quantile([0.01,0.99]))