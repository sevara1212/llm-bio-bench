import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')
adata_full = adata.raw.to_adata()
adata_full.obs['leiden'] = adata.obs['leiden']

markers = ['CD3D','CD3E','IL7R','CCR7','CD4','CD8A','CD8B','MS4A1','CD79A','LYZ','CD14','FCGR3A','MS4A7',
           'NKG7','GNLY','FCER1A','CST3','PPBP','PF4','GZMA','GZMB','CCL5']

X = adata_full[:, [g for g in markers if g in adata_full.var_names]].X
X = np.asarray(X.todense()) if hasattr(X, 'todense') else np.asarray(X)
present = [g for g in markers if g in adata_full.var_names]
df = pd.DataFrame(X, columns=present, index=adata_full.obs_names)
df['leiden'] = adata_full.obs['leiden'].values
means = df.groupby('leiden').mean()
pd.set_option('display.width', 250)
print(means.round(2))