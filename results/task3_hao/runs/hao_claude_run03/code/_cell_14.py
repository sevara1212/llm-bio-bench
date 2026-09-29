import scanpy as sc
import pandas as pd

full = sc.read_h5ad('full_with_markers.h5ad')
genes = ['CD1C','FCER1A','CLEC10A','ITGAX','SIRPA','CD68','MRC1']
genes_present=[g for g in genes if g in full.var_names]
print(genes_present)
df = sc.get.obs_df(full, keys=genes_present+['leiden'])
mean_expr = df.groupby('leiden').mean()
print(mean_expr.loc[['6','7','13','17']].round(2).T)