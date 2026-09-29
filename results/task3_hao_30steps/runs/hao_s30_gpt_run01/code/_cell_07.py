import scanpy as sc
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.6'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=50)
for c in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=c)
 genes=[x for x in d.names if not x.startswith(('RPL','RPS','MT-'))][:15]
 print(f'{c:>2} {sum(ad.obs[key]==c):4}:', ', '.join(genes))