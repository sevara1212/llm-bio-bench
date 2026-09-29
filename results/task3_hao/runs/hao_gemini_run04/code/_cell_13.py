import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')

# Check expression of key lineage markers across all clusters
markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'CCR7', 'CD14', 'FCGR3A', 'MS4A1', 'CD79A', 'GNLY', 'NKG7', 'NCAM1', 'CST3', 'CLEC4C', 'LILRA4', 'PPBP', 'PF4', 'FCER1A', 'CLEC9A', 'CD1C', 'MZB1', 'JCHAIN', 'STMN1', 'MKI67', 'CDK6', 'GZMK', 'FOXP3']

df = sc.get.obs_df(adata, keys=['leiden_0.5'] + [m for m in markers if m in adata.raw.var_names])
mean_expr = df.groupby('leiden_0.5').mean()
print("Mean expression of markers across clusters:")
pd.set_option('display.max_columns', 35)
pd.set_option('display.width', 1000)
print(mean_expr)