import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('full_with_clusters.h5ad')

markers = ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','SELL','TCF7',
           'NKG7','GNLY','KLRD1','KLRF1','NCAM1',
           'MS4A1','CD19','CD79A','CD79B','TCL1A','IGHD','IGHM',
           'MZB1','JCHAIN','TNFRSF17',
           'CD14','LYZ','FCN1','S100A8','S100A9','VCAN',
           'FCGR3A','MS4A7','LST1',
           'CST3','HLA-DRA','CLEC9A','CD1C','FCER1A',
           'IL3RA','LILRA4','CLEC4C','GZMB',
           'PPBP','PF4','ITGA2B',
           'HBB','HBA1',
           'MKI67','STMN1','TOP2A',
           'PRSS57','CD34']
markers = [m for m in markers if m in adata.var_names]

df = sc.get.obs_df(adata, keys=markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
mean_expr.to_csv('mean_expr_by_cluster.csv')
print("saved")
print(mean_expr.shape)