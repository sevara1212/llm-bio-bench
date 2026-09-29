import scanpy as sc,pandas as pd
ad=sc.read_h5ad('processed.h5ad')
# use log norm raw (not scaled)
sc.tl.rank_genes_groups(ad,'leiden_0.4',method='wilcoxon',use_raw=True,n_genes=15)
r=ad.uns['rank_genes_groups']
for c in ad.obs.leiden_0.4.cat.categories:
 print('\nCLUSTER',c,'n=',sum(ad.obs.leiden_0.4==c))
 print(', '.join(r['names'][c][:15]))