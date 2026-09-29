import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()

# Basic QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)

n_before = adata.n_obs
# Filter cells
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.pct_counts_mt < 15].copy()
n_after_cells = adata.n_obs

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)
n_genes_after = adata.n_vars

print(f'cells before: {n_before}, after: {n_after_cells}')
print(f'genes after filter: {n_genes_after}')

adata.write('qc_filtered.h5ad')