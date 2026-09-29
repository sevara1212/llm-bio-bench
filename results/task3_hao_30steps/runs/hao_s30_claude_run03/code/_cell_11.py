import scanpy as sc
import pandas as pd
import numpy as np

ad = sc.read_h5ad('clustered.h5ad')

markers = ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','SELL',
    'GZMK','GZMB','GZMH','NKG7','GNLY','KLRD1','KLRF1','FCGR3A',
    'MS4A1','CD79A','CD79B','IGHD','TCL1A','MZB1','JCHAIN','TNFRSF17',
    'CD14','LYZ','S100A8','S100A9','FCN1',
    'MS4A7','LST1','AIF1',
    'CST3','FCER1A','CD1C','CLEC9A','BATF3',
    'LILRA4','IL3RA','IRF7','TCF4',
    'PPBP','PF4','NRGN',
    'HBB','HBA1',
    'STMN1','MKI67','PRSS57','CD34']

df = sc.get.obs_df(ad, keys=markers+['leiden_05'], use_raw=True)
mean_expr = df.groupby('leiden_05').mean()
pd.set_option('display.width', 300)
pd.set_option('display.max_columns', 60)

# print in chunks of 10 columns
cols = mean_expr.columns.tolist()
for i in range(0, len(cols), 8):
    print(mean_expr[cols[i:i+8]].round(2))
    print()