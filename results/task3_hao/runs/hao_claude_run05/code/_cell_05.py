import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("Before filtering:", adata.shape)

# Filter cells
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs['pct_counts_mt'] < 15].copy()
adata = adata[adata.obs['n_genes_by_counts'] < 6000].copy()

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)

print("After filtering:", adata.shape)
adata.write('adata_filtered.h5ad')