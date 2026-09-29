import scanpy as sc
import numpy as np
adata = sc.read_h5ad('raw_counts.h5ad')

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], inplace=True, percent_top=None)
print(adata.obs[['n_genes_by_counts','total_counts','pct_counts_mt']].describe())