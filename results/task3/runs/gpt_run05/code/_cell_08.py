import scanpy as sc
A=sc.read_h5ad('pbmc_processed_prelim.h5ad'); A.obs['cluster']=A.obs['leiden_1.0']
sc.tl.rank_genes_groups(A,'cluster',method='wilcoxon',use_raw=False,n_genes=30)
d=sc.get.rank_genes_groups_df(A,group=None)
for g in sorted(d.group.unique(),key=int):
 print('\n',g,'N',(A.obs.cluster==g).sum(),':',', '.join(d[d.group==g].head(25).names.tolist()))