import scanpy as sc,pandas as pd,numpy as np
ad=sc.read_h5ad('pbmc_processed_hvg.h5ad'); ad.obs['cluster']=ad.obs['leiden_0.8'].astype(str)
print(ad.obs.groupby('cluster')[['n_genes_by_counts','total_counts','pct_counts_mt']].median().round(1))
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=True,n_genes=500)
r=sc.get.rank_genes_groups_df(ad,group=None)
for c in ['0','6']:
 x=r[r.group==c];x=x[~x.names.str.match(r'^(RP[SL]|MT-|MALAT1|EEF|LDHB|TPT1|B2M)')]
 print('\n',c,x[['names','scores']].head(40).to_string(index=False))