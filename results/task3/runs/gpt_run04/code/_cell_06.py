import scanpy as sc, pandas as pd
ad=sc.read_h5ad('pbmc_processed.h5ad')
for key in ['leiden_0.5','leiden_0.7','leiden_1.0']:
 sc.tl.rank_genes_groups(ad,groupby=key,method='wilcoxon',use_raw=False,pts=True,key_added='r')
 print('\n###',key)
 d=sc.get.rank_genes_groups_df(ad,group=None,key='r')
 for c in ad.obs[key].cat.categories:
  x=d[d.group==c].head(12)
  print(c, ', '.join(x.names.tolist()))