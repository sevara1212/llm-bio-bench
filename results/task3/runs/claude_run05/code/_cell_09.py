import scanpy as sc
adata_full = sc.read_h5ad('adata_clustered.h5ad').raw.to_adata()
adata = sc.read_h5ad('adata_clustered.h5ad')
adata_full.obs['leiden'] = adata.obs['leiden']

# Check canonical markers across clusters
markers = ['CD3D','CD3E','CD4','CD8A','IL7R','CCR7','S100A4','MS4A1','CD79A','CD14','LYZ','FCGR3A','MS4A7',
           'GNLY','NKG7','GZMB','FCER1A','CST3','PPBP','PF4']
import pandas as pd
df = sc.get.obs_df(adata_full, keys=['leiden']+markers)
means = df.groupby('leiden').mean()
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
print(means.round(2))