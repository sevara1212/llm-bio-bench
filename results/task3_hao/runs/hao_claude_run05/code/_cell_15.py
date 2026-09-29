import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_final.h5ad')

# double check FOXP3 in cluster 16
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden'].values
for g in ['FOXP3','CD8A','IL2RA']:
    if g in adata_raw.var_names:
        vals = pd.Series(adata_raw[:, g].X.toarray().flatten(), index=adata_raw.obs_names)
        vals = vals.groupby(adata_raw.obs['leiden']).mean()
        print(g)
        print(vals.round(2))