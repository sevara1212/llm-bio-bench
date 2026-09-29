import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Check if standard QC filtering should be applied
# Let's see: min_genes=200, min_cells=3, pct_counts_mt < 15 or 20
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print("Before QC:", adata.shape)

# Let's see how many cells are filtered if pct_counts_mt < 15 vs no cell filtering
# Note min_genes in dataset is 499, min total_counts is 898, max mt is 15.06%
# If someone prepared raw_counts.h5ad for a benchmark, they usually either:
# 1) expect you to do standard QC: sc.pp.filter_cells(adata, min_genes=200), sc.pp.filter_genes(adata, min_cells=3)
# But wait, what if they expect labels for all 7841 barcodes?
# Let's check if any cell had < 200 genes: 0 cells!
# If we filter genes with min_cells=3:
sc.pp.filter_genes(adata, min_cells=3)
print("After gene filtering (min_cells=3):", adata.shape)
# All 7841 cells are preserved!