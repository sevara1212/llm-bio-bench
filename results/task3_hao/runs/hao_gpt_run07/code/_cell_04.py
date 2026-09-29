import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_counts.h5ad'); ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-'); sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True)
for g in [500,5500,6000]: print(g,(ad.obs.n_genes_by_counts>=g).sum(),(ad.obs.n_genes_by_counts>g).sum())
for t in [40000,50000,60000]: print(t,(ad.obs.total_counts>t).sum())
print(ad.obs.sort_values('n_genes_by_counts',ascending=False)[['total_counts','n_genes_by_counts','pct_counts_mt']].head(20))