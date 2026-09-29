import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_counts.h5ad');ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-');sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True)
keep=(ad.obs.n_genes_by_counts>=500)&(ad.obs.n_genes_by_counts<=5500)&(ad.obs.pct_counts_mt<12)
print(keep.sum(), 'removed',(~keep).sum())
print('reasons low/high/mt', (ad.obs.n_genes_by_counts<500).sum(),(ad.obs.n_genes_by_counts>5500).sum(),(ad.obs.pct_counts_mt>=12).sum())
print(ad.obs.loc[~keep,['total_counts','n_genes_by_counts','pct_counts_mt']].describe())