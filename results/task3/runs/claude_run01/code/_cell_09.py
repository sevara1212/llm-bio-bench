import scanpy as sc
import numpy as np
adata = sc.read_h5ad('clustered.h5ad')
adata_full = adata.raw.to_adata()
adata_full.obs['leiden'] = adata.obs['leiden']

markers = ['CD3D','CD3E','IL7R','CCR7','CD4','CD8A','CD8B','MS4A1','CD79A','LYZ','CD14','FCGR3A','MS4A7',
           'NKG7','GNLY','FCER1A','CST3','PPBP','PF4','GZMA','GZMB','CCL5']

import pandas as pd
df = pd.DataFrame(index=adata_full.obs['leiden'].cat.categories)
for g in markers:
    if g in adata_full.var_names:
        vals = adata_full[:, g].X
        vals = np.asarray(vals).flatten()
        means = pd.Series(vals, index=adata_full.obs_names).groupby(adata_full.obs['leiden']).mean()
        df[g] = means
pd.set_option('display.width', 200)
print(df.round(2))