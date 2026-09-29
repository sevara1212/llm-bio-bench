import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad')
sc.tl.rank_genes_groups(ad,'leiden_0.6',method='wilcoxon',use_raw=True,n_genes=12)
for c in ad.obs['leiden_0.6'].cat.categories:
 d=sc.get.rank_genes_groups_df(ad,group=c).head(10)
 print('\nCLUSTER',c,'n=',sum(ad.obs['leiden_0.6']==c))
 print(', '.join(d.names.tolist()))
# canonical mean expression all clusters
markers=['CD3D','CD3E','TRAC','IL7R','LTB','CCR7','MALAT1','NKG7','GNLY','PRF1','GZMB','KLRD1','MS4A1','CD79A','CD74','HLA-DRA','CD37','CD14','LYZ','S100A8','S100A9','FCGR3A','LST1','CTSS','FCER1A','CST3','CLEC10A','GZMB','JCHAIN','MZB1','SDC1','PPBP','PF4','ITGA2B','FCGR3B','CSF3R','TPSAB1','KIT']
# raw logged mean
for c in ad.obs['leiden_0.6'].cat.categories:
 x=ad[ad.obs['leiden_0.6']==c].raw.to_adata()[:,[g for g in markers if g in ad.var_names]].X.mean(axis=0)
 import numpy as np
 x=np.asarray(x).ravel(); gg=[g for g in markers if g in ad.var_names]
 tops=sorted(zip(x,gg),reverse=True)[:8]
 print(c, ' '.join(f'{g}:{v:.1f}' for v,g in tops))