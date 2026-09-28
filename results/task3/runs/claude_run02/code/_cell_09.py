import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_ranked.h5ad')

extra_markers = ['CCR7','SELL','S100A4','CD27','TCF7','GZMK','CD3D','IL7R','LTB','AQP3']
extra_markers = [m for m in extra_markers if m in adata.var_names]

df = sc.get.obs_df(adata, keys=extra_markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
print(mean_expr.round(2))

print(adata.obs['leiden'].value_counts())