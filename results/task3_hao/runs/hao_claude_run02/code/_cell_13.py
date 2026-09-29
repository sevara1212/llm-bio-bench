import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('clustered2.h5ad')

markers = ['GZMB','GZMH','NKG7','GNLY','KLRD1','KLRF1','FCGR3A',
           'CD14','LYZ','S100A8','S100A9','FCN1','VCAN','MS4A7','CD68',
           'MS4A1','CD79A','CD79B','CD19','IGHD','IGHM','TCL1A','BANK1',
           'MZB1','JCHAIN','TNFRSF17',
           'PPBP','PF4','GNG11','TUBB1',
           'HBB','HBA1','ALAS2',
           'IL3RA','LILRA4','CLEC9A','BATF3','CD1C','FCER1A']
markers = [m for m in markers if m in adata.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden'], use_raw=False)
mean_expr = df.groupby('leiden').mean().round(2)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
cols = mean_expr.columns.tolist()
for i in range(0, len(cols), 8):
    print(mean_expr[cols[i:i+8]])
    print()