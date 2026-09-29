import scanpy as sc
ad=sc.read_h5ad('processed_markers.h5ad');key='leiden_0.6'
for cl in ad.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=cl)
 # filter generic
 bad=('RPL','RPS','MT-','MALAT1','HBA','HBB','B2M','TMSB','EEF','GAPDH','ACTB')
 genes=[x for x in d.names if not x.startswith(bad)][:20]
 print(f'{cl:>2} {(ad.obs[key]==cl).sum():>4}: '+', '.join(genes))