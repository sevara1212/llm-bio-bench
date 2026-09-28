import scanpy as sc
A=sc.read_h5ad('pbmc_processed_prelim.h5ad')
for key in ['leiden_0.4','leiden_0.5','leiden_0.6','leiden_0.7','leiden_0.8','leiden_1.0']:
 print(key,A.obs[key].value_counts().sort_index().to_dict())
# markers at .8
A.obs['cluster']=A.obs['leiden_0.8']
sc.tl.rank_genes_groups(A,'cluster',method='wilcoxon',use_raw=False,n_genes=20)
d=sc.get.rank_genes_groups_df(A,group=None)
for g in sorted(d.group.unique(),key=int):
 print('\n',g,'N',(A.obs.cluster==g).sum(),':',', '.join(d[d.group==g].head(18).names.tolist()))