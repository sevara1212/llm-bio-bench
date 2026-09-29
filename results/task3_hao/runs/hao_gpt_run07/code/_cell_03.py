import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
ad.var['ribo']=ad.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True,percent_top=[20])
for c in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 print(c, np.percentile(ad.obs[c],[0,1,5,25,50,75,95,99,100]))
print('mt genes',ad.var.mt.sum())
print('cells min genes', (ad.obs.n_genes_by_counts<200).sum(), 'mt >20', (ad.obs.pct_counts_mt>20).sum())