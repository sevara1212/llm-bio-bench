import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('pbmc_processed.h5ad'); ad.obs['cluster']=ad.obs['leiden_0.8'].astype(str)
# derive expression from raw counts retained layer
ad.X=ad.layers['counts'].copy(); sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',n_genes=25)
r=sc.get.rank_genes_groups_df(ad,group=None)
for c in map(str,range(8)):
 print(c, (ad.obs.cluster==c).sum(), ':', ', '.join(r[r.group==c].names.head(20).tolist()))
# inspect canonical exact mean log expr
G=['IL7R','CCR7','LTB','MAL','TCF7','LEF1','CD3D','TRAC','LST1','LYZ','S100A8','S100A9','CTSS','FCGR3A','MS4A7','LGALS3','IFITM3','MS4A1','CD79A','CD37','CD74','HLA-DRA','NKG7','GNLY','PRF1','GZMB','CCL5','FCER1A','CD1C','CLEC10A','PPBP','PF4']
means=pd.DataFrame(index=sorted(ad.obs.cluster.unique(),key=int))
for g in G:
 means[g]=[float(ad[ad.obs.cluster==c,g].X.mean()) if g in ad.var_names else 0 for c in means.index]
print(means.round(2).to_string())
# save markers
r.to_csv('cluster_markers.csv',index=False)
ad.write('pbmc_clustered.h5ad')