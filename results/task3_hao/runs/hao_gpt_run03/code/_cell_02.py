import scanpy as sc, os
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head()); print(ad.var.head())
print(ad.X.dtype, ad.X.shape, 'sparse', hasattr(ad.X,'tocsr'))
print(ad.obs.columns.tolist(),ad.var.columns.tolist())
print('layers',list(ad.layers), 'raw',ad.raw is not None)