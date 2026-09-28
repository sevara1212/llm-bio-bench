import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# standard QC filters for this dataset (classic scanpy tutorial thresholds)
print("Before filter:", adata.shape)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()

print("After filter:", adata.shape)

adata.write('qc.h5ad')