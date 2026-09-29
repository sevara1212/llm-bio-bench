import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')

# Check mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
print("Number of MT genes:", adata.var['mt'].sum())
# Check ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
print("Number of Ribo genes:", adata.var['ribo'].sum())

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
print(adata.obs.describe())