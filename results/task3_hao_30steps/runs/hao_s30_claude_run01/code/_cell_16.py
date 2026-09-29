import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

genes_check = ['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2',
               'CD4','CD8A','CD8B','FOXP3','IL2RA','CTLA4',
               'GNLY','NKG7','KLRD1','KLRB1','KLRF1','NCAM1','FCGR3A','GZMK','GZMB','GZMH','GZMA','CCL5','CST7','PRF1',
               'SELL','CCR7','TCF7','LEF1','IL7R','CD27','CD28']

genes_present = [g for g in genes_check if g in adata.var_names]
expr = adata[:, genes_present].X
if not isinstance(expr, np.ndarray):
    expr = expr.toarray()
df2 = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df2['leiden'] = adata.obs['leiden'].values

means = df2.groupby('leiden').mean()
print(means.loc[['0','1','2','4','5','8','16','18','20']].round(2).T)