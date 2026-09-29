import scanpy as sc, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
for sep in ['_','-']:
 print(sep, ad.obs_names.str.split(sep).str[0].value_counts())