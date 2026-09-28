import scanpy as sc, os, numpy as np, pandas as pd
p='filtered_gene_bc_matrices/hg19'
a=sc.read_10x_mtx(p, var_names='gene_symbols', cache=False)
print(a)
print(a.var.head(15)); print(a.obs.head()); print('dups',a.var_names.duplicated().sum())
print('total range',np.asarray(a.X.sum(1)).ravel().min(),np.percentile(np.asarray(a.X.sum(1)).ravel(),[1,5,50,95,99]),np.asarray(a.X.sum(1)).ravel().max())
print('genes detected',np.percentile(np.diff(a.X.indptr),[0,1,5,50,95,99,100]))