import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True, log1p=False)
print(ad.obs[['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']].describe(percentiles=[.01,.05,.1,.5,.9,.95,.99]).round(2))
print('mito genes',ad.var.mt.sum(),'ribo',ad.var.ribo.sum())
print('zeros genes', (ad.var.n_cells_by_counts==0).sum())
# quantile high and preliminary filters outcomes
for g in [100,200,300,500]:
 for mt in [10,15,20,25]: print(g,mt,((ad.obs.n_genes_by_counts>=g)&(ad.obs.pct_counts_mt<mt)).sum())