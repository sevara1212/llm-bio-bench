import scanpy as sc
import pandas as pd
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)

ad = sc.read_h5ad('processed.h5ad')

genes = ['CD3D','CD3E','TRAC','CD4','CD8A','CD8B','GZMK','CD27','TCF7','KLRB1','CD1C','FCER1A','CLEC9A','BATF3','XCR1']
genes = [g for g in genes if g in ad.raw.var_names]
df = sc.get.obs_df(ad, keys=genes+['leiden'], use_raw=True)
means = df.groupby('leiden').mean()
print(means.round(2).to_string())