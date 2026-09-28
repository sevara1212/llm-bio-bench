import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("Before filtering:", adata.shape)

# Basic filtering: genes expressed in >=3 cells, cells with >=200 genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print("After basic filter:", adata.shape)

# Filter on QC metrics (standard PBMC3k thresholds)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

print("After QC filter:", adata.shape)

adata.write('qc_filtered.h5ad')