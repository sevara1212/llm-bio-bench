import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('adata_full_lognorm_clustered.h5ad')

markers = ['CD8A','CD8B','CD4','CD3D','GZMK','SLC4A10','KLRB1','ZBTB16','TRAV1-2','CCR7','SELL','TCF7','LEF1','FOXP3','IL2RA','CTLA4','IKZF2']
markers = [m for m in markers if m in adata.var_names]
df = sc.get.obs_df(adata, keys=markers+['leiden'])
mean_expr = df.groupby('leiden').mean()
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 60)
print(mean_expr.round(3).T)