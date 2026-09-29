import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

genes_check = ['SLC4A10','ZBTB16','KLRB1','TRDV1','TRDV2','TRGV9','NCAM1','FCGR3A','KLRC1','XCL1','XCL2','SPON2','FGFBP2']

genes_present = [g for g in genes_check if g in adata.var_names]
expr = adata[:, genes_present].X
if not isinstance(expr, np.ndarray):
    expr = expr.toarray()
df2 = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df2['leiden'] = adata.obs['leiden'].values

means = df2.groupby('leiden').mean()
print(means.loc[['4','8','16','20']].round(2).T)