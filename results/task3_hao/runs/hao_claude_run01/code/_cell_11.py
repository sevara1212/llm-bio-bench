import scanpy as sc
import pandas as pd
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)

ad = sc.read_h5ad('processed.h5ad')

for name, genes in {
    'Tnaive_mem': ['IL7R','CCR7','SELL','TCF7'],
    'NK': ['NKG7','GNLY','KLRD1','KLRB1','GZMK','GZMB','FCGR3A'],
}.items():
    genes = [g for g in genes if g in ad.raw.var_names]
    df = sc.get.obs_df(ad, keys=genes+['leiden'], use_raw=True)
    means = df.groupby('leiden').mean()
    print(f"--- {name} ---")
    print(means.round(2).to_string())
    print()