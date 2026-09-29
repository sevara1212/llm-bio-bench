import scanpy as sc
import numpy as np
ad = sc.read_h5ad('raw_counts.h5ad')

ad.var['mt'] = ad.var_names.str.startswith('MT-')
ad.var['ribo'] = ad.var_names.str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
print(ad.obs[['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo']].describe())