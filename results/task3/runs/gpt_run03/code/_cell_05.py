import scanpy as sc, pandas as pd
A=sc.read_h5ad('pbmc_processed.h5ad')
for r in [.7,.9,1.1]:
 sc.tl.leiden(A,resolution=r,key_added=f'l{r}',random_state=0)
 print(r,A.obs[f'l{r}'].value_counts().sort_index().to_dict())
# resolution .9 markers
A.obs['leiden']=A.obs['l0.9'].astype(str)
sc.tl.rank_genes_groups(A,'leiden',method='wilcoxon',use_raw=False,n_genes=40)
for cl in sorted(A.obs.leiden.unique(),key=int):
 d=sc.get.rank_genes_groups_df(A,group=cl).head(18)
 print('\n',cl,len(A.obs[A.obs.leiden==cl]));print(' '.join(d.names))
A.write('pbmc_processed_r09.h5ad')