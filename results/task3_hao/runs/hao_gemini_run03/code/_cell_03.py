import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

adata = sc.read_h5ad('raw_counts.h5ad')
print(type(adata.X))
print("X shape:", adata.X.shape)
if sp.issparse(adata.X):
    print("X min, max:", adata.X.data.min(), adata.X.data.max())
else:
    print("X min, max:", np.min(adata.X), np.max(adata.X))

# Check mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
print("Number of MT genes:", adata.var['mt'].sum())
# Check ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
print("Number of Ribo genes:", adata.var['ribo'].sum())

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())