import scanpy as sc
adata=sc.read_h5ad('pbmc_qc_clustered.h5ad'); key='leiden_0.5'
# rankings not saved prior. recompute
sc.tl.rank_genes_groups(adata,key,method='wilcoxon',use_raw=False,n_genes=15)
for cl in adata.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(adata,group=cl,key='rank_genes_groups')
 print(f'{cl}\t{(adata.obs[key]==cl).sum()}\t'+' '.join(d.names.head(10)))