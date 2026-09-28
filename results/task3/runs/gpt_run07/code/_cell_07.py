import scanpy as sc,pandas as pd,numpy as np
ad=sc.read_h5ad('pbmc_processed_hvg.h5ad');ad.obs['cluster']=ad.obs['leiden_0.8'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=True,n_genes=100)
r=sc.get.rank_genes_groups_df(ad,group=None)
for c in sorted(ad.obs.cluster.unique(),key=int):
 print('\n',c,ad.obs.cluster.value_counts()[c], ', '.join(r[r.group==c].names.head(30).tolist()))
markers=['IL7R','CCR7','LTB','MALAT1','LST1','S100A8','S100A9','LYZ','FCGR3A','MS4A7','NKG7','GNLY','PRF1','GZMB','TRBC2','CD3D','CD3E','IL32','MS4A1','CD79A','CD37','CD74','HLA-DRA','CD14','PPBP','PF4','FCER1A','CST3']
D={}
for g in markers:
 if g in ad.raw.var_names:
  x=ad.raw[:,g].X.toarray().ravel()
  D[g]=[x[(ad.obs.cluster==c).to_numpy()].mean() for c in sorted(ad.obs.cluster.unique(),key=int)]
print(pd.DataFrame(D,index=sorted(ad.obs.cluster.unique(),key=int)).round(2).to_string())