import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp

adata = sc.read_h5ad('adata_processed.h5ad')

markers = ['CD3D','CD3E','TRAC','IL7R','CD4','CD8A','CD8B','CCR7','GZMK','GZMB','NKG7','GNLY','KLRD1','KLRF1',
           'MS4A1','CD79A','CD79B','BANK1','CD14','LYZ','S100A8','S100A9','FCN1','FCGR3A','MS4A7',
           'CST3','FCER1A','CLEC10A','CLEC9A','LILRA4','IL3RA','CLEC4C','PPBP','PF4','GP9',
           'HBB','HBA1','MZB1','JCHAIN','MKI67','STMN1','TYMS','FOXP3','IL2RA']

genes_present = [g for g in markers if g in adata.raw.var_names]
expr = adata.raw[:, genes_present].X
if sp.issparse(expr):
    expr = expr.toarray()
df = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df['leiden'] = adata.obs['leiden'].values
mean_expr = df.groupby('leiden').mean()
mean_expr.T.to_csv('marker_mean_expr.csv')

# print in chunks
cols = mean_expr.T.columns.tolist()
for i in range(0,len(cols),6):
    print(mean_expr.T.iloc[:, i:i+6].round(2))
    print()