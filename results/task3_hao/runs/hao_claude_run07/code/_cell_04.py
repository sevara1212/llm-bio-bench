import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

import numpy as np
for q in [0.01,0.02,0.05,0.95,0.98,0.99]:
    print(q, np.quantile(adata.obs['n_genes_by_counts'],q), np.quantile(adata.obs['total_counts'],q), np.quantile(adata.obs['pct_counts_mt'],q))

print((adata.obs['pct_counts_mt']>10).sum())
print((adata.obs['n_genes_by_counts']<200).sum())