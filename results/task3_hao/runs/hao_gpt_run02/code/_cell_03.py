import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True,log1p=False)
for c in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 print('\n',c, ad.obs[c].quantile([0,.01,.05,.1,.25,.5,.75,.9,.95,.99,1]).to_string())
print('mito genes',ad.var.mt.sum())
print('cells gene cut', [(x,int((ad.obs.n_genes_by_counts>=x).sum())) for x in [100,200,300,500,800,1000]])
print('mt cut',[(x,int((ad.obs.pct_counts_mt<x).sum())) for x in [10,15,20,25,30]])