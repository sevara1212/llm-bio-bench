import scanpy as sc
ad = sc.read_h5ad('raw_counts.h5ad')
print(ad)
print(ad.X[:10,:10])
print(ad.obs.head())
print(ad.var.head())