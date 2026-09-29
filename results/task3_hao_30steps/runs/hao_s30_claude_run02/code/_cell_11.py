import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')

markers = ['NKG7','GNLY','KLRD1','KLRF1','NCAM1','FCGR3A',
           'MS4A1','CD79A','CD79B','CD19','IGHD','TCL1A','MZB1','JCHAIN',
           'CD14','LYZ','S100A8','S100A9','FCN1','VCAN','MS4A7','CST3']
markers = [m for m in markers if m in adata.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 60)
print(mean_expr.round(2).T)