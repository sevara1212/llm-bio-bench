import scanpy as sc, pandas as pd
A=sc.read_h5ad('pbmc_processed_prelim.h5ad')
print(A.obs.leiden.value_counts().sort_index())
df=sc.get.rank_genes_groups_df(A,group=None)
for g in sorted(df.group.unique(),key=int):
 x=df[df.group==g].head(20)
 print('\nCLUSTER',g,'n=',(A.obs.leiden==g).sum())
 print(', '.join(x.names.tolist()))