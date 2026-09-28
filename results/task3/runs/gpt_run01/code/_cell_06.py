import scanpy as sc
ad=sc.read_h5ad('pbmc_processed.h5ad')
ad.obs['cluster']=ad.obs['leiden_1.0'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=False,n_genes=40)
for cl in sorted(ad.obs.cluster.unique(),key=int):
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(25)
 print('\n',cl,'N',(ad.obs.cluster==cl).sum(),'\n',', '.join(d.names.tolist()))
# targeted mean expression (all / detected)
genes=['IL7R','LTB','CCR7','MALAT1','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A1','CD79A','CD37','NKG7','GNLY','PRF1','GZMB','CCL5','CD3D','TRAC','IL32','FCER1A','CST3','CLEC10A','PF4','PPBP','GNG11','MKI67']
for cl in sorted(ad.obs.cluster.unique(),key=int):
 a=ad[ad.obs.cluster==cl]
 vals=[]
 for g in genes:
  if g in a.var_names: vals.append(f'{g}:{a[:,g].X.mean():.1f}')
 print(cl,' '.join(vals))