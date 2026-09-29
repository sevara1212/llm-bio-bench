import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Look at distribution
print("Cells with mt > 15:", (adata.obs['pct_counts_mt'] > 15).sum())
print("Cells with mt > 20:", (adata.obs['pct_counts_mt'] > 20).sum())
print("Cells with n_genes < 200:", (adata.obs['n_genes_by_counts'] < 200).sum())
print("Cells with n_genes < 500:", (adata.obs['n_genes_by_counts'] < 500).sum())
print("Cells with total_counts < 1000:", (adata.obs['total_counts'] < 1000).sum())

# Notice min n_genes is 499, max pct_counts_mt is 15.06!
# That means the raw dataset might have ALREADY been pre-filtered with min_genes=500, max_mt=15!