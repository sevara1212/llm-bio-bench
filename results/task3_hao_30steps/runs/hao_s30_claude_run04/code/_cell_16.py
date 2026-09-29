import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_full_annotated.h5ad')

genes_check2 = ['GZMK','GZMB','GZMH','NKG7','GNLY','KLRD1','FCGR3A','NCAM1',
               'MS4A1','CD79A','CD79B','CD19','IGHD','TCL1A',
               'MZB1','JCHAIN','TNFRSF17',
               'CD14','LYZ','S100A8','FCN1','VCAN',
               'MS4A7','LST1',
               'CST3','FCER1A','CLEC9A','CD1C']

df = sc.get.obs_df(adata, keys=genes_check2+['leiden_r05'])
mean_expr = df.groupby('leiden_r05').mean()
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 250)
print(mean_expr.round(2).T)