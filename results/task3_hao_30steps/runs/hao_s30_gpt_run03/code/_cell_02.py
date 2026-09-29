import scanpy as sc, os
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head()); print(ad.var.head()); print(ad.layers.keys()); print(ad.X.__class__)
print(ad.obs.columns.tolist(), ad.var.columns.tolist())
print('sum',ad.X.sum())