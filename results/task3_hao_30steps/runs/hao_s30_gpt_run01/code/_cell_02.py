import scanpy as sc, os
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head()); print(ad.var.head())
print(ad.X.dtype, ad.X.min(),ad.X.max())
print(ad.obs.columns.tolist(), ad.var.columns.tolist())
print('sum', ad.X.sum(axis=1)[:5])