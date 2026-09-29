import scanpy as sc, os
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head()); print(ad.var.head())
print(ad.X.dtype, ad.X.shape, 'nnz', ad.X.nnz if hasattr(ad.X,'nnz') else 'dense')
print(ad.obs.columns.tolist(), ad.var.columns.tolist())
print('total', ad.X.sum())