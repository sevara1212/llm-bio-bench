import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True, log1p=False)
for x in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 print('\n',x, ad.obs[x].describe(percentiles=[.01,.05,.1,.25,.5,.75,.9,.95,.99]).round(2))
print('mt genes',sum(ad.var.mt))
print('threshold remain',((ad.obs.n_genes_by_counts>=200)&(ad.obs.n_genes_by_counts<=6000)&(ad.obs.pct_counts_mt<15)).sum())