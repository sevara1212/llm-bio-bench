import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')

extra = ['CD4','CD8A','CD8B','TRDC','TRGC1','TRGC2','TRDV1','TRDV2','CD3D','FOXP3','IL2RA','CTLA4','TIGIT','GZMK','CD27','TCF7','CCR7','SELL','KLRG1','GNLY']
extra = [g for g in extra if g in adata.raw.var_names]
df = sc.get.obs_df(adata, keys=extra+['leiden_05'], use_raw=True)
m = df.groupby('leiden_05').mean()
pd.set_option('display.width', 300)
pd.set_option('display.max_columns', 30)
print(m.loc[['17','18']].round(3))