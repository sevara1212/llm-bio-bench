import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad')
sc.tl.rank_genes_groups(ad,'leiden_0.6',method='wilcoxon',use_raw=True,pts=True)
for cl in ad.obs.leiden_0.6.cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(12)
 print('\nCLUSTER',cl,'n=',(ad.obs.leiden_0.6==cl).sum())
 print(', '.join(d.names.tolist()))
ad.write('processed_markers.h5ad')