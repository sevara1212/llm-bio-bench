import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# QC filtering thresholds
print("Before filtering:", adata.shape)

# filter cells
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]

# filter genes expressed in at least 3 cells
sc.pp.filter_genes(adata, min_cells=3)

print("After filtering:", adata.shape)

adata.write('qc_filtered.h5ad')