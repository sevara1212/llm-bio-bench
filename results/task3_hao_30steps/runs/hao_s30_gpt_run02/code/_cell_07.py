import scanpy as sc
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.6'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,pts=True)
for cl in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(12)
 print('\nCLUSTER',cl,'n=',(ad.obs[key]==cl).sum())
 print(', '.join(d.names.tolist()))
ad.write('processed_markers.h5ad')