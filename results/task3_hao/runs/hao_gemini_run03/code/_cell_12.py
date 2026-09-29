import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')
markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'LEF1', 'IL7R', 'S100A4', 
           'FOXP3', 'IL2RA', 'TRDC', 'NKG7', 'GNLY', 'FCGR3A',
           'MS4A1', 'CD19', 'CD79A', 'IGHM', 'IGHD', 'CD27', 'MZB1',
           'CD14', 'LYZ', 'S100A9', 'MS4A7', 'FCER1A', 'CST3', 'CLEC9A', 'IL3RA', 'LILRA4',
           'PPBP', 'PF4', 'HBB', 'HBA1', 'MKI67', 'CD34']

avail = [m for m in markers if m in adata.var_names]
df_expr = pd.DataFrame(index=[str(i) for i in range(21)])
for m in avail:
    df_expr[m] = [float(adata[adata.obs['leiden_0.5'] == str(i), m].X.mean()) for i in range(21)]

pd.set_option('display.max_columns', 40)
pd.set_option('display.width', 1000)
print(df_expr.round(2))