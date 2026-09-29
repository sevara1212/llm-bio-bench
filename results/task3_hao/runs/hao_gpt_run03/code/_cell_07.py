import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad')
# markers use log normalized raw slot (not regressed/scaled)
sc.tl.rank_genes_groups(ad,'leiden_0.5',method='wilcoxon',use_raw=True,n_genes=15)
res=ad.uns['rank_genes_groups']
for c in ad.obs.leiden_0.5.cat.categories:
 d=pd.DataFrame({k:res[k][c] for k in ['names','scores','pvals_adj']}).head(10)
 print('\nCLUSTER',c,'n=',(ad.obs.leiden_0.5==c).sum()); print(', '.join(d.names.astype(str)))
ad.write_h5ad('processed_markers.h5ad')