import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
print(adata)
print("Any duplicated var_names?", adata.var_names.duplicated().sum())
adata.var_names_make_unique()
print(adata)