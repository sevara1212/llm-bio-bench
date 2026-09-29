import scanpy as sc,numpy as np,pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad'); ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-'); ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True,log1p=False)
for x in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 q=np.percentile(ad.obs[x],[0,1,5,10,25,50,75,90,95,99,100]);print(x,np.round(q,2))
print('mt genes',ad.var.mt.sum())
# top cell count values
print(ad.obs[['total_counts','n_genes_by_counts','pct_counts_mt']].describe())