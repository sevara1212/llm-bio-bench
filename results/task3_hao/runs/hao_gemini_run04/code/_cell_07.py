import pandas as pd
adata = sc.read_h5ad('raw_counts.h5ad')
print("Total cells in raw:", len(adata.obs_names))
print("Any duplicates?", adata.obs_names.duplicated().sum())