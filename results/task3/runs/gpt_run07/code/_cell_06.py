import scanpy as sc,pandas as pd
ad=sc.read_h5ad('pbmc_processed_hvg.h5ad')
ad.obs['cluster']=ad.obs['leiden_0.8'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=True,n_genes=20)
r=sc.get.rank_genes_groups_df(ad,group=None)
for c in sorted(ad.obs.cluster.unique(),key=int):
 print('\nCLUSTER',c,'n=',sum(ad.obs.cluster==c))
 print(r[r.group==c][['names','scores','pvals_adj']].head(15).to_string(index=False))
# canonical marker means counts from log raw
markers=['IL7R','LTB','CCR7','MALAT1','LST1','S100A8','S100A9','LYZ','FCGR3A','MS4A7','NKG7','GNLY','PRF1','GZMB','TRBC2','CD3D','CD3E','CCR7','LTB','IL32','MS4A1','CD79A','CD37','CD74','HLA-DRA','CD14','PPBP','PF4','FCER1A','CST3']
D=pd.DataFrame(index=sorted(ad.obs.cluster.unique(),key=int))
for g in markers:
 if g in ad.raw.var_names:
  x=ad.raw[:,g].X
  D[g]=[x[ad.obs.cluster==c].mean() for c in D.index]
print('\nmeans');print(D.round(2).to_string())