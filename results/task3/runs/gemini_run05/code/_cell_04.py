import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=True)
print(adata)
print(adata.obs.head())
print(adata.var.head())