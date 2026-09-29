import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.5'; sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=100)
r=ad.uns['rank_genes_groups']
for c in ad.obs[key].cat.categories:
 genes=[str(x) for x in r['names'][c] if not str(x).startswith(('RPL','RPS','MT-'))][:18]
 print(c, (ad.obs[key]==c).sum(), ':', ', '.join(genes))