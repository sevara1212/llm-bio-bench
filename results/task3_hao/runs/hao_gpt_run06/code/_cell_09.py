import scanpy as sc
adata=sc.read_h5ad('pbmc_qc_clustered.h5ad');key='leiden_0.5';sc.tl.rank_genes_groups(adata,key,method='wilcoxon',use_raw=False,n_genes=40)
for cl in ['0','4','6','8','9','10','13','15','16','17','19','20']:
 d=sc.get.rank_genes_groups_df(adata,group=cl,key='rank_genes_groups')
 print('\n',cl, ', '.join(d.names.head(35)))