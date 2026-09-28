import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('pbmc_processed.h5ad'); key='leiden_1.0'
# rank on raw log-normalized? process X scaled, layer unavailable normalized. Can assess group means on scaled markers.
markers=['CD3D','CD3E','TRAC','IL7R','LTB','CCR7','MALAT1','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A1','CD79A','CD74','NKG7','GNLY','GZMB','GZMK','CCL5','FCER1A','CST3','CLEC10A','PPBP','PF4','MKI67','TYMS','TOP2A','CD14','HLA-DRA']
for c in ad.obs[key].cat.categories:
 x=ad[ad.obs[key]==c,:]
 means=np.asarray(x[:,markers].X.mean(axis=0)).ravel()
 print('\n',c,'n=',x.n_obs)
 print(' '.join(f'{g}:{v:.1f}' for g,v in zip(markers,means) if v>.25))
# Proper differential stored use saved recomputation
sc.tl.rank_genes_groups(ad,groupby=key,method='wilcoxon',key_added='rank')
d=sc.get.rank_genes_groups_df(ad,None,key='rank')
for c in ad.obs[key].cat.categories:
 print(c, d[d.group==c].head(30).names.tolist())