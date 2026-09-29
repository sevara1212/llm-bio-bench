import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')

markers = ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','SELL','TCF7',
           'CD14','LYZ','FCGR3A','MS4A7','CST3',
           'MS4A1','CD79A','CD79B','IGHD','IGHM','TCL1A',
           'NKG7','GNLY','KLRD1','GZMB','GZMK',
           'PPBP','PF4','ITGA2B',
           'FCER1A','CLEC10A','CLEC9A','IL3RA','LILRA4',
           'MKI67','STMN1',
           'MZB1','JCHAIN',
           'HBB','HBA1',
           'PRSS57','CD34','KLRB1','CD27','FOXP3','IKZF2']

markers = [m for m in markers if m in adata.raw.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden_05'], use_raw=True)
mean_expr = df.groupby('leiden_05').mean()
mean_expr.round(2).to_csv('mean_expr.csv')
print("saved")
print(mean_expr.iloc[:, 12:24].round(2))