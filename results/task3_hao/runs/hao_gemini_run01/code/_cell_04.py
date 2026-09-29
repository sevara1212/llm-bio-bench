import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('raw_counts.h5ad')

# Check mitochondrial genes
mito_genes = adata.var_names.str.startswith('MT-')
print(f"Number of MT genes: {mito_genes.sum()}")
print("MT genes:", adata.var_names[mito_genes].tolist()[:10])

# Ribosomal genes
ribo_genes = adata.var_names.str.startswith(('RPS', 'RPL'))
print(f"Number of ribosomal genes: {ribo_genes.sum()}")

adata.var['mt'] = mito_genes
adata.var['ribo'] = ribo_genes

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

print(adata.obs.describe().T[['min', '25%', '50%', '75%', 'max']])