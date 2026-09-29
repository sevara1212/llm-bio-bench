import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('clustered2.h5ad')

markers = ['MS4A1','CD79A','CD79B','IGHD','IGHM','TCL1A','BANK1','MZB1','JCHAIN','TNFRSF17','PPBP','PF4','GNG11','TUBB1','LYZ','S100A8','S100A9','FCN1','VCAN','MS4A7','CD68']
markers = [m for m in markers if m in adata.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden'], use_raw=False)
mean_expr = df.groupby('leiden').mean().round(2)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
cols = mean_expr.columns.tolist()
for i in range(0, len(cols), 7):
    print(mean_expr[cols[i:i+7]])
    print()