import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=False)
print("Adata shape:", adata.shape)
print("Var names sample:", adata.var_names[:10])
print("Obs names sample:", adata.obs_names[:10])