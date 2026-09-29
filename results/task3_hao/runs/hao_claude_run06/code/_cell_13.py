import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('clustered_final.h5ad')

markers = ['CD3D','CD3E','CD4','CD8A','CD8B','CCR7','SELL','IL7R','TCF7',
           'NKG7','GNLY','GZMB','KLRD1','NCAM1','FCGR3A',
           'MS4A1','CD79A','CD79B','MZB1','JCHAIN',
           'LYZ','CD14','FCN1','MS4A7','CST3',
           'LILRA4','IL3RA','PPBP','PF4','HBB','HBA1',
           'MKI67','STMN1','PRSS57','CD34','CLEC9A','FCER1A']

markers = [m for m in markers if m in adata.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden_final'])
mean_expr = df.groupby('leiden_final').mean()
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 50)
print(mean_expr.round(2))