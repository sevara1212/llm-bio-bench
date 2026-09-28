import scanpy as sc
adata_full = sc.read_h5ad('processed_full.h5ad')
adata_hvg = sc.read_h5ad('clustered.h5ad')
adata_full.obs['leiden'] = adata_hvg.obs['leiden'].values

extra = ['SELL','CCR7','TCF7','LEF1','GPR183','FOXP3','IL2RA','KLRB1']
extra = [g for g in extra if g in adata_full.var_names]
import pandas as pd
pd.set_option('display.max_columns', None)
df = sc.get.obs_df(adata_full, keys=extra+['leiden'])
print(df.groupby('leiden').mean().round(2))