import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed_pbmc.h5ad'); key='leiden_0.6'
# use raw expression in .raw
x=ad.raw.to_adata(); x.obs[key]=ad.obs[key]
# test ranks based original expression
sc.tl.rank_genes_groups(x,key,method='wilcoxon',n_genes=50)
for c in x.obs[key].cat.categories:
 print('\nCLUSTER',c,'N=',sum(x.obs[key]==c))
 print(', '.join(x.uns['rank_genes_groups']['names'][c][:30]))
# lineage expression means
mk=['CD3D','CD4','CD8A','CD8B','CCR7','LTB','IL7R','MALAT1','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A7','NKG7','KLRD1','GNLY','GZMB','CCL5','MS4A1','CD79A','CD74','HLA-DRA','TCL1A','MZB1','JCHAIN','PF4','PPBP','GATA3','TRDC','MKI67','TYMS']
df=pd.DataFrame(x[:,mk].X.toarray(),index=x.obs[key]).groupby(level=0,observed=True).mean()
print(df.round(2).to_string())