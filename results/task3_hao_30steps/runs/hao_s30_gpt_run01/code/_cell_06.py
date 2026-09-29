import scanpy as sc
ad=sc.read_h5ad('processed.h5ad')
sc.tl.rank_genes_groups(ad,'leiden_0.6',method='wilcoxon',use_raw=True,n_genes=30)
for c in ad.obs.leiden_0.6.cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=c)
 # exclude ribos/mt
 genes=[x for x in d.names if not x.startswith(('RPL','RPS','MT-'))][:15]
 print(f'{c:>2} {sum(ad.obs.leiden_0.6==c):4}:', ', '.join(genes))