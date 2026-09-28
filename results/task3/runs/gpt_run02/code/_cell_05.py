import scanpy as sc, pandas as pd
ad=sc.read_h5ad('pbmc_processed.h5ad'); ad.obs['cluster']=ad.obs['leiden_0.8'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=False,n_genes=30)
r=sc.get.rank_genes_groups_df(ad,group=None)
for c in sorted(ad.obs.cluster.unique(),key=int):
 print('\nCLUSTER',c,'n=',(ad.obs.cluster==c).sum())
 print(r[r.group==c][['names','scores','logfoldchanges','pvals_adj']].head(15).to_string(index=False))
# marker score means
markers={'B':['MS4A1','CD79A','CD37','HLA-DRA','CD74','CD79B','CD83'],'CD4T':['IL7R','LTB','CCR7','MALAT1','LST1','MAL','LEF1','TCF7','CCR7'],'CD8T_NK':['NKG7','GNLY','GZMB','PRF1','CTSW','CD3D','TRAC','CCL5'],'Mono':['LYZ','S100A8','S100A9','LGALS3','CTSD','FCN1','TYROBP','LST1'],'DC':['FCER1A','CST3','CD1C','CLEC10A'],'Platelet':['PPBP','PF4','NRGN','GNG11']}
for k,gs in markers.items():
 present=[x for x in gs if x in ad.var_names]
 x=ad[:,present].X.mean(axis=1).A1 if hasattr(ad[:,present].X.mean(axis=1),'A1') else ad[:,present].X.mean(axis=1)
 print(k, pd.Series(x,index=ad.obs.cluster).groupby(level=0).mean().round(2).to_dict())