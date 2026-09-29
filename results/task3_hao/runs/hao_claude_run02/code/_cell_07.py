import scanpy as sc
adata_hvg = sc.read_h5ad('clustered.h5ad')  # has neighbors etc but let's redo from adata_hvg saved separately better