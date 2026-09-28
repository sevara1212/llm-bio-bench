import scanpy as sc
ad=sc.read_h5ad('pbmc_processed.h5ad');ad.obs['cluster']=ad.obs['leiden_1.0'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=False,n_genes=100)
for cl in ['1','7']:
 print('\n###',cl)
 d=sc.get.rank_genes_groups_df(ad,group=cl)
 print(d.head(100).to_string(index=False))