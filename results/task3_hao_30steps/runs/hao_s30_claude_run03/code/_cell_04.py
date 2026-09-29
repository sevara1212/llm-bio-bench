import scanpy as sc
import numpy as np

ad = sc.read_h5ad('raw_counts.h5ad')
ad.var['mt'] = ad.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

for q in [0.01,0.02,0.05,0.95,0.98,0.99]:
    print(q, 'genes', np.quantile(ad.obs.n_genes_by_counts,q), 'counts', np.quantile(ad.obs.total_counts,q), 'mt', np.quantile(ad.obs.pct_counts_mt,q))