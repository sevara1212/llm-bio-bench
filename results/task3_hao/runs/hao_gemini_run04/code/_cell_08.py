import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
print(f"Initial shape: {adata.shape}")

# Quality control
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Keep track of raw barcodes
all_barcodes = adata.obs_names.copy()

# Filter cells: min_genes=200, pct_counts_mt < 15% (almost all cells pass, min is 499)
# Let's filter out mt > 15% (only 3 cells) or should we keep all cells in labels.csv?
# Let's check how many cells are filtered if we do min_genes=200, pct_mt < 15%
sc.pp.filter_genes(adata, min_cells=3)
print(f"After gene filtering: {adata.shape}")