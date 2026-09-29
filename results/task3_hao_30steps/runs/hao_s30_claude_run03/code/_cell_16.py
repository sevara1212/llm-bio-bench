import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('clustered.h5ad')
extra_markers = ['FOXP3','IL2RA','CD8A','CD8B','CD4','CCR7','TCF7','SELL','SELPLG','KLRB1','CD27','GPR183']
df = sc.get.obs_df(ad, keys=extra_markers+['leiden_05'], use_raw=True)
mean_expr = df.groupby('leiden_05').mean()
print(mean_expr.round(2))