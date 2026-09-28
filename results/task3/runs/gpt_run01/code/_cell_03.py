import scanpy as sc, pandas as pd, numpy as np
p='filtered_gene_bc_matrices/hg19'
adata=sc.read_10x_mtx(p, var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
print(adata)
print(adata.var.head())
print(adata.obs.head())
print('total',adata.X.sum(), 'range', np.asarray(adata.X.sum(1)).ravel().min(),np.asarray(adata.X.sum(1)).ravel().max())
print('genes',adata.var_names[:20].tolist())