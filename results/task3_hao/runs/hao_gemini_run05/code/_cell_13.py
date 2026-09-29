import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('processed.h5ad')

# Key marker expression per cluster
genes_to_check = [
    'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'IL7R', 'S100A4', 'FOXP3',
    'NKG7', 'GNLY', 'NCAM1', 'FCGR3A',
    'MS4A1', 'CD19', 'CD79A', 'MZB1', 'JCHAIN',
    'CD14', 'LYZ', 'VCAN', 'MS4A7',
    'FCER1A', 'CST3', 'CLEC9A', 'CLEC10A', 'IL3RA', 'LILRA4', 'TCF4',
    'PPBP', 'PF4',
    'MKI67', 'STMN1', 'TOP2A'
]

genes_in_data = [g for g in genes_to_check if g in adata.var_names]
df_expr = pd.DataFrame(index=adata.obs['leiden'].cat.categories)

for g in genes_in_data:
    df_expr[g] = [adata[adata.obs['leiden'] == c, g].raw.X.mean() for c in df_expr.index]

pd.set_option('display.max_columns', 35)
pd.set_option('display.width', 1000)
print(df_expr.round(2))