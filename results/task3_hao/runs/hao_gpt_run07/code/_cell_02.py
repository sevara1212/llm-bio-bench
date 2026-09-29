import scanpy as sc, os
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head()); print(ad.var.head())
print(ad.X.shape, ad.X.dtype)
print('layers',ad.layers.keys()); print('raw',ad.raw)
print('obs cols',ad.obs.columns.tolist())