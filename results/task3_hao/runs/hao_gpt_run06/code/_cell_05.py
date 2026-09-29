import scanpy as sc, pandas as pd
adata=sc.read_h5ad('pbmc_qc_clustered.h5ad')
key='leiden_0.5'
sc.tl.rank_genes_groups(adata,key,method='wilcoxon',use_raw=False,n_genes=20,pts=True)
for cl in adata.obs[key].cat.categories:
 d=sc.get.rank_genes_groups_df(adata,group=cl, key='rank_genes_groups')
 print('\nCLUSTER',cl,'n=',(adata.obs[key]==cl).sum())
 print(', '.join(d.names.head(12).tolist()))