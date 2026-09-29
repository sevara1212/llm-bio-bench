import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("n_genes_by_counts quantiles:")
print(adata.obs['n_genes_by_counts'].quantile([0.01,0.02,0.05,0.95,0.98,0.99,0.995,1.0]))
print("total_counts quantiles:")
print(adata.obs['total_counts'].quantile([0.01,0.02,0.05,0.95,0.98,0.99,0.995,1.0]))
print("pct_counts_mt quantiles:")
print(adata.obs['pct_counts_mt'].quantile([0.5,0.9,0.95,0.98,0.99,1.0]))

# genes detected in how many cells
print((adata.var['n_cells_by_counts']==0).sum() if 'n_cells_by_counts' in adata.var else "no col")