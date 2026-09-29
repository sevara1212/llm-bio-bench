import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_full_annotated.h5ad')

genes_check = ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','SELL','TCF7',
               'GZMK','GZMB','GZMH','NKG7','GNLY','KLRD1','FCGR3A','NCAM1',
               'MS4A1','CD79A','CD79B','CD19','IGHD','TCL1A',
               'MZB1','JCHAIN','TNFRSF17',
               'CD14','LYZ','S100A8','FCN1','VCAN',
               'MS4A7','LST1',
               'CST3','FCER1A','CLEC9A','CD1C',
               'IL3RA','IRF7','IRF8',
               'PPBP','PF4','GNG11',
               'HBB','HBA1','HBA2',
               'FOXP3','IL2RA','TNFRSF18','TNFRSF4',
               'STMN1','MKI67',
               'PRSS57','CD34','SOX4']

df = sc.get.obs_df(adata, keys=genes_check+['leiden_r05'])
mean_expr = df.groupby('leiden_r05').mean()
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 250)
print(mean_expr.round(2).T)