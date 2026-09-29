import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)

print(np.histogram(adata.obs.pct_counts_mt, bins=20))
print('n_genes_by_counts<200', (adata.obs.n_genes_by_counts<200).sum())
print('total_counts<500', (adata.obs.total_counts<500).sum())