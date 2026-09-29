import scanpy as sc
ad=sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.obs.head(),ad.obs.columns.tolist())
print(ad.var.head(),ad.var.columns.tolist())
print('X',ad.X.min(),ad.X.max(),ad.X.sum())