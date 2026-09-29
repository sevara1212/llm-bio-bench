import scanpy as sc
import numpy as np

ad = sc.read_h5ad('raw_counts.h5ad')
ad.var['mt'] = ad.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print("Before filtering:", ad.shape)

# Basic QC filtering
sc.pp.filter_cells(ad, min_genes=200)
sc.pp.filter_genes(ad, min_cells=3)

# filter high mito / outlier cells
ad = ad[ad.obs.pct_counts_mt < 15, :].copy()
ad = ad[ad.obs.n_genes_by_counts < 5500, :].copy()  # remove potential doublets/outliers
ad = ad[ad.obs.total_counts < 40000, :].copy()

print("After filtering:", ad.shape)
ad.write('qc_filtered.h5ad')