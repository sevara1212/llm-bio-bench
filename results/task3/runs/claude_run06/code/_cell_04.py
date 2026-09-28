import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw.h5ad')

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(adata.obs[['n_genes_by_counts','total_counts','pct_counts_mt']].describe())