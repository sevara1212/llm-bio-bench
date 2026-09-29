import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# check quantiles for upper bounds
for q in [0.9,0.95,0.99,0.995,0.999]:
    print(q, adata.obs['n_genes_by_counts'].quantile(q), adata.obs['total_counts'].quantile(q), adata.obs['pct_counts_mt'].quantile(q))

print("n cells before filter", adata.n_obs)