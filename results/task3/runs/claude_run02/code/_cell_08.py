import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('adata_ranked.h5ad')

markers = ['CD3D','CD3E','IL7R','CD4','CD8A','CD8B','CCR7','SELL',
           'CD14','LYZ','FCGR3A','MS4A7','FCN1',
           'MS4A1','CD79A','CD79B',
           'NKG7','GNLY','GZMB','NCAM1',
           'FCER1A','CLEC10A','CST3',
           'PPBP','PF4','GP9']

markers = [m for m in markers if m in adata.var_names]

df = sc.get.obs_df(adata, keys=markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 30)
print(mean_expr.round(2))