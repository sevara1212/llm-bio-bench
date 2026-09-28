import scanpy as sc
import pandas as pd
import numpy as np

# Load 10x dataset
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
print(adata)
print("Barcodes head:", adata.obs_names[:5])
print("Genes head:", adata.var_names[:5])