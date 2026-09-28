import scanpy as sc
adata = sc.read_h5ad('clustered.h5ad')
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden']

genes_check = ['CD3D','CD3E','IL7R','CD4','CD8A','CD8B','CCR7','S100A4',
               'CD14','LYZ','FCGR3A','MS4A7','FCER1A','CLEC10A',
               'MS4A1','CD79A','CD79B',
               'GNLY','NKG7','KLRB1','GZMB','GZMK','CCL5',
               'PPBP','PF4']
import pandas as pd
df = sc.get.obs_df(adata_raw, keys=genes_check+['leiden'])
print(df.groupby('leiden').mean().T)