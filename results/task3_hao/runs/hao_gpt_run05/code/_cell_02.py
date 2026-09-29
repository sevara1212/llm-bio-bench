import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head());print(ad.var.head())
print('layers',ad.layers.keys(),'raw',ad.raw)
print(ad.X.max(),ad.X.min(), type(ad.X))
print('genes',ad.var_names[:20].tolist())