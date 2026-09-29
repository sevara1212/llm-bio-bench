import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

genes_check = ['IGHD','IGHM','TCL1A','CD27','CD38','AIM2','TNFRSF13B','SELL','CD24']
genes_present = [g for g in genes_check if g in adata.var_names]
expr = adata[:, genes_present].X
if not isinstance(expr, np.ndarray):
    expr = expr.toarray()
df2 = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df2['leiden'] = adata.obs['leiden'].values
means = df2.groupby('leiden').mean()
print(means.loc[['3','10']].T)