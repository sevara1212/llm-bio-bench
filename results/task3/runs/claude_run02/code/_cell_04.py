import matplotlib
matplotlib.use('Agg')
import scanpy as sc
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# filter genes present in few cells
sc.pp.filter_genes(adata, min_cells=3)

# Standard PBMC3k-like QC thresholds
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
adata = adata[adata.obs.n_genes_by_counts > 200, :]

print(adata.shape)

adata.write('adata_filtered.h5ad')