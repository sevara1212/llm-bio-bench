import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
for name,pat in [('mt','^MT-'),('ribo','^RP[SL]'),('hb','^HB[AB]')]: ad.var[name]=ad.var_names.str.match(pat)
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo','hb'],inplace=True,percent_top=[50,100,200])
print(ad.obs[['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo','pct_counts_hb']].describe(percentiles=[.01,.05,.1,.25,.5,.75,.9,.95,.99]))
print('genes mt',ad.var.mt.sum(),'ribo',ad.var.ribo.sum(),'hb',ad.var.hb.sum())
print('zero genes', (ad.var.n_cells_by_counts==0).sum())
# Top genes
print(ad.var.sort_values('total_counts',ascending=False)[['total_counts','n_cells_by_counts','mt','ribo','hb']].head(25))