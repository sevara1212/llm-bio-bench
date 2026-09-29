import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# quantiles for potential upper thresholds
for q in [0.01,0.02,0.5,0.9,0.95,0.98,0.99,0.995,0.999]:
    print(q, adata.obs['total_counts'].quantile(q), adata.obs['n_genes_by_counts'].quantile(q))
print('mt pct quantiles')
for q in [0.9,0.95,0.98,0.99]:
    print(q, adata.obs['pct_counts_mt'].quantile(q))