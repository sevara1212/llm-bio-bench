import scanpy as sc
import pandas as pd

full = sc.read_h5ad('full_with_markers.h5ad')
genes = ['CD14','FCGR3A','MS4A7','FCN1','S100A8','S100A9','LYZ','CST3','CD74','CLEC9A','BATF3','LILRA4','IL3RA','IRF7']
df = sc.get.obs_df(full, keys=genes+['leiden'])
mean_expr = df.groupby('leiden').mean()
print(mean_expr.loc[['6','7','13','17','9']].round(2).T)