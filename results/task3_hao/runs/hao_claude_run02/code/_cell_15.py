import scanpy as sc
import pandas as pd
adata = sc.read_h5ad('clustered2.h5ad')

genes = ['CD3D','CD3E','CD3G','TRDC','TRGC1','TRGC2','KLRB1','IL7R','TCF7','SELL','CD4','CD8A','CD8B','FOXP3','CTLA4','IKZF2','NKG7','GNLY','KLRD1']
genes = [g for g in genes if g in adata.var_names]
df = sc.get.obs_df(adata, keys=genes+['leiden'], use_raw=False)

# fraction expressing (>0)
frac = df.groupby('leiden').apply(lambda x: (x[genes]>0).mean())
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
print(frac.round(2).loc[['5','6','7','14','15','16','18','19']])