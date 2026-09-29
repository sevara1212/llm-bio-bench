import scanpy as sc
import pandas as pd
import numpy as np

# Load 10x dataset
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
print(adata)
print("Shape:", adata.shape)
print("First 5 barcodes:", adata.obs_names[:5].tolist())
print("First 5 genes:", adata.var_names[:5].tolist())