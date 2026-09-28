import scanpy as sc, numpy as np
ad=sc.read_h5ad('pbmc_processed.h5ad'); key='leiden_1.0'
markers=[x for x in ['CD3D','CD3E','IL7R','LTB','CCR7','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A1','CD79A','CD74','NKG7','GNLY','GZMB','GZMK','CCL5','FCER1A','CST3','PPBP','PF4','MKI67','TYMS','TOP2A','CD14','HLA-DRA'] if x in ad.var_names]
for c in ad.obs[key].cat.categories:
 x=ad[ad.obs[key]==c,:]
 means=np.asarray(x[:,markers].X.mean(axis=0)).ravel()
 print(c, 'n=',x.n_obs, '|', ' '.join(f'{g}:{v:.1f}' for g,v in zip(markers,means) if v>.25))
sc.tl.rank_genes_groups(ad,groupby=key,method='wilcoxon',key_added='rank')
d=sc.get.rank_genes_groups_df(ad,None,key='rank')
for c in ad.obs[key].cat.categories:
 print('\n',c, d[d.group==c].head(25).names.tolist())