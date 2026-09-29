import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True,log1p=False)
for c in ['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_ribo']:
 print('\n',c, ad.obs[c].quantile([0,.01,.05,.1,.25,.5,.75,.9,.95,.99,1]).round(2).to_dict())
print('mt genes',ad.var.mt.sum(), 'zero genes',(ad.var.total_counts==0).sum())
print('cells <200',sum(ad.obs.n_genes_by_counts<200),'mt>20',sum(ad.obs.pct_counts_mt>20),'mt>15',sum(ad.obs.pct_counts_mt>15),'genes>6000',sum(ad.obs.n_genes_by_counts>6000))