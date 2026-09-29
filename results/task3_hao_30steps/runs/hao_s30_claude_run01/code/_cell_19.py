import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('full_with_clusters.h5ad')

genes_check = ['FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF18','TNFRSF4','CD3D','CD3E','TRAC','TRDC','KLRB1','IL7R','SELL','CCR7']
genes_present = [g for g in genes_check if g in adata.var_names]
expr = adata[:, genes_present].X
if not isinstance(expr, np.ndarray):
    expr = expr.toarray()
df2 = pd.DataFrame(expr, columns=genes_present, index=adata.obs_names)
df2['leiden'] = adata.obs['leiden'].values

means = df2.groupby('leiden').mean()
print(means.loc[['20']].T)

print()
print('n cells cluster20:', (adata.obs.leiden=='20').sum())
print('n_genes_by_counts cluster20 mean:', adata.obs.loc[adata.obs.leiden=='20','n_genes_by_counts'].mean())
print('total_counts cluster20 mean:', adata.obs.loc[adata.obs.leiden=='20','total_counts'].mean())
print('overall mean n_genes:', adata.obs['n_genes_by_counts'].mean())