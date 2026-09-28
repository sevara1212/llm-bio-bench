import scanpy as sc, pandas as pd, numpy as np
p='filtered_gene_bc_matrices/hg19'
adata=sc.read_10x_mtx(p, var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
print(adata)
print(adata.var.head()); print(adata.obs.head())
print('counts:', adata.X.sum(), 'cells',adata.n_obs,'genes',adata.n_vars)
print('gene names examples', adata.var_names[:20].tolist())