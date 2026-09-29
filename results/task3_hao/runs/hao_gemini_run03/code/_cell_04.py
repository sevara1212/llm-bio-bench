import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']].describe())

# Check percentiles
for col in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']:
    print(col, np.percentile(adata.obs[col], [0.1, 1, 5, 50, 95, 99, 99.5, 99.9, 100]))