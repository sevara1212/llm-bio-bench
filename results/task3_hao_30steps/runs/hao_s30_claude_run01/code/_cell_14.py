import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

genes_check = ['CD3D','CD3E','CD4','CD8A','CD8B','FOXP3','IL2RA','TRDC','TRGC1','TRGC2',
               'MKI67','TOP2A','CD34','AVP','KIT','GATA2',
               'AXL','SIGLEC6','CLEC9A','CD1C','FCER1A','ITGAX','BATF3']

genes_present = [g for g in genes_check if g in adata.var_names]
expr = adata[:, genes_present].X
if not isinstance(expr, np.ndarray):
    expr = expr.toarray()
df2 = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df2['leiden'] = adata.obs['leiden'].values

means = df2.groupby('leiden').mean()
print(means.loc[['4','12','15','18','20','23']].round(2).T)