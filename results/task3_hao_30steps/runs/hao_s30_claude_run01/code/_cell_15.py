import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

genes_check = ['CD14','FCGR3A','LYZ','S100A8','S100A9','VCAN','FCN1','MS4A7','LST1','CDKN1C',
               'CLEC9A','BATF3','CD1C','FCER1A','CST3','HLA-DRA','ITGAX',
               'NLRP3','IL1B','NFKBIA','G0S2']

genes_present = [g for g in genes_check if g in adata.var_names]
expr = adata[:, genes_present].X
if not isinstance(expr, np.ndarray):
    expr = expr.toarray()
df2 = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df2['leiden'] = adata.obs['leiden'].values

means = df2.groupby('leiden').mean()
print(means.loc[['6','7','17','19','21']].round(2).T)