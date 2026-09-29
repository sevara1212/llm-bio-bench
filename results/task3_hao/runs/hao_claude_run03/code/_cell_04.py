import scanpy as sc
import numpy as np

adata = sc.read_h5ad('raw_counts.h5ad')
adata.var_names_make_unique()
print(adata.obs_names[:5].tolist())

# check for batch/prefix info in barcodes
prefixes = set([b.split('_')[0] for b in adata.obs_names])
print(prefixes)