import scanpy as sc
import pandas as pd

raw_adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
raw_adata.var_names_make_unique()
print("Raw barcodes count:", raw_adata.n_obs)

processed_adata = sc.read_h5ad('processed_adata.h5ad')
print("Processed barcodes count:", processed_adata.n_obs)