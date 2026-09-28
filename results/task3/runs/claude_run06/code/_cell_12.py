import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')
raw = sc.read_h5ad('raw.h5ad')

print("Raw cells:", raw.n_obs)
print("QC-passed cells:", adata.n_obs)
print("Filtered out:", raw.n_obs - adata.n_obs)