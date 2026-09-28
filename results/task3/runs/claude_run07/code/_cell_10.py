import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)

import scanpy as sc
adata_full = sc.read_h5ad('processed_full.h5ad')
adata_hvg = sc.read_h5ad('clustered.h5ad')
adata_full.obs['leiden'] = adata_hvg.obs['leiden'].values

markers = ['CD3D','CD3E','CD4','IL7R','CCR7','S100A4','CD8A','CD8B','GZMK','NKG7','GNLY',
           'MS4A1','CD79A','CD14','LYZ','FCGR3A','MS4A7','FCER1A','CST3','PPBP','PF4']
markers = [m for m in markers if m in adata_full.var_names]

df = sc.get.obs_df(adata_full, keys=markers+['leiden'])
means = df.groupby('leiden').mean().round(2)
print(means)

print()
print(adata_full.obs['leiden'].value_counts())