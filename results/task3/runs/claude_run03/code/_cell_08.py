import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('adata_clustered.h5ad')
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden']

# Check canonical marker genes across clusters
markers = ['CD3D','CD3E','IL7R','CD4','CD8A','CD8B','CCR7','SELL',
           'CD14','LYZ','FCN1','FCGR3A','MS4A7',
           'CD79A','MS4A1','CD79B',
           'GNLY','NKG7','GZMB','KLRD1',
           'FCER1A','CD1C','CLEC10A',
           'PPBP','PF4','GP9']

df = sc.get.obs_df(adata_raw, keys=markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 30)
print(mean_expr.round(2))