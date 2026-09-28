import scanpy as sc
ad=sc.read_h5ad('pbmc_processed.h5ad');ad.obs['cluster']=ad.obs.leiden_1.0.astype(str);sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',n_genes=60)
for cl in ['1','7']:
 d=sc.get.rank_genes_groups_df(ad,group=cl)
 print(cl, list(d.names[:60]))