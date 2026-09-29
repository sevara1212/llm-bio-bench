import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.5'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=15)
res=ad.uns['rank_genes_groups']
for c in ad.obs[key].cat.categories:
 d=pd.DataFrame({k:res[k][c] for k in ['names','scores','pvals_adj']}).head(10)
 print('\nCLUSTER',c,'n=',(ad.obs[key]==c).sum()); print(', '.join(d.names.astype(str)))
ad.write_h5ad('processed_markers.h5ad')