import scanpy as sc, numpy as np
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.6'
for c in ['4','5','7','12','19']:
 sub=ad[ad.obs[key]==c]
 # use rank broader
 sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=200)
 d=sc.get.rank_genes_groups_df(ad,group=c)
 g=[x for x in d.names if not x.startswith(('RPL','RPS','MT-'))][:45]
 print('\n',c,', '.join(g))