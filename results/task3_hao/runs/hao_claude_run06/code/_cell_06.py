import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

n0 = adata.shape[0]
print("Before filtering:", adata.shape)

# Filter genes expressed in fewer than 3 cells
sc.pp.filter_genes(adata, min_cells=3)
print("After gene filter:", adata.shape)

# Filter cells: standard PBMC thresholds
adata = adata[(adata.obs['n_genes_by_counts'] >= 200) &
              (adata.obs['n_genes_by_counts'] <= 6000) &
              (adata.obs['total_counts'] >= 500) &
              (adata.obs['pct_counts_mt'] < 15)].copy()
print("After cell filter:", adata.shape)
print("Cells removed:", n0 - adata.shape[0])

adata.write('filtered.h5ad')