import scanpy as sc,pandas as pd
ad=sc.read_h5ad('pbmc_processed.h5ad');ad.obs['cluster']=ad.obs['leiden_1.0'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=False,n_genes=200)
for cl in ['0','1','7']:
 d=sc.get.rank_genes_groups_df(ad,group=cl)
 print('\n',cl)
 print(d[(d['scores']>5)&(~d.names.str.startswith(('RPS','RPL','MT-')))].head(60).to_string(index=False))