import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

n_before = adata.n_obs
# QC filters
adata = adata[adata.obs['n_genes_by_counts'] > 200, :]
adata = adata[adata.obs['n_genes_by_counts'] < 6000, :]
adata = adata[adata.obs['pct_counts_mt'] < 15, :]
adata = adata[adata.obs['total_counts'] < 40000, :]

sc.pp.filter_genes(adata, min_cells=3)

print('before', n_before, 'after', adata.n_obs, 'genes', adata.n_vars)
adata.write('qc_filtered.h5ad')