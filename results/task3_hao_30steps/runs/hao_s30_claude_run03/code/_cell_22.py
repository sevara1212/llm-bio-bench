import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('final_labeled.h5ad')

# Sanity check - look at raw barcode index and confirm it matches original raw data
raw = sc.read_h5ad('raw_counts.h5ad')
print("Filtered barcodes:", ad.n_obs, "of raw:", raw.n_obs)
print(ad.obs_names[:5].tolist())
print(raw.obs_names[:5].tolist())

# check that this is a subset
print(set(ad.obs_names).issubset(set(raw.obs_names)))