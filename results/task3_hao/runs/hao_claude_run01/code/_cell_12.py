import scanpy as sc
import pandas as pd
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)

ad = sc.read_h5ad('processed.h5ad')

genes = ['CD3D','CD4','CD8A','FOXP3','IL2RA','CTLA4','TNFRSF18','TNFRSF4','KLF2','SELL','CCR7',
         'CD34','PRSS57','SOX4','MPO','GATA1',
         'PPBP','HBB',
         'LILRA4','GZMB','IRF8','CLEC9A']
genes = [g for g in genes if g in ad.raw.var_names]
df = sc.get.obs_df(ad, keys=genes+['leiden'], use_raw=True)
means = df.groupby('leiden').mean()
print(means.round(2).to_string())

print(ad.obs['leiden'].value_counts().sort_index())