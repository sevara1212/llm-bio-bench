import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# check percentiles for upper bounds
for col in ['n_genes_by_counts','total_counts','pct_counts_mt']:
    print(col, np.percentile(adata.obs[col], [1,5,95,99,99.5,100]))