import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True,log1p=False)
print(ad.obs[['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']].describe(percentiles=[.01,.05,.1,.25,.5,.75,.9,.95,.99]))
print('MT genes',ad.var.mt.sum())
print('top counts'); print(ad.obs.sort_values('total_counts',ascending=False)[['total_counts','n_genes_by_counts','pct_counts_mt']].head(15))