import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
for key, pat in [('mt','MT-'),('ribo','^RP[SL]'),('hb','^HB[AB]')]: ad.var[key]=ad.var_names.str.startswith(pat) if key=='mt' else ad.var_names.str.contains(pat,regex=True)
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo','hb'],inplace=True,log1p=True)
print(ad.obs[['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo','pct_counts_hb']].describe(percentiles=[.01,.05,.1,.25,.5,.75,.9,.95,.99]).round(2))
print('mt genes',ad.var.mt.sum())
print('top genes', ad.var.sort_values('n_cells_by_counts',ascending=False)[['n_cells_by_counts','mean_counts']].head(20))
ad.obs.to_csv('qc_metrics.csv')