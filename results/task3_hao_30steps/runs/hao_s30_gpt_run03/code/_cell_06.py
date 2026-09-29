import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad')
key='leiden_0.6'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=15)
r=ad.uns['rank_genes_groups']
for cl in ad.obs[key].cat.categories:
 print('\n',cl,ad.obs[key].value_counts()[cl],':', ', '.join(r['names'][cl][:12]))
# canonical marker means
markers=['CD3D','CD3E','TRAC','IL7R','LTB','CCR7','MALAT1','NKG7','GNLY','PRF1','GZMB','KLRD1','MS4A1','CD79A','CD74','HLA-DRA','CD14','LYZ','S100A8','S100A9','FCGR3A','LST1','FCER1A','CLEC10A','GZMB','JCHAIN','MZB1','MKI67','PPBP','PF4','GATA2','TPSAB1','HBB']
df=sc.get.obs_df(ad,keys=[key]+[x for x in markers if x in ad.var_names]).groupby(key,observed=True).mean().T
pd.set_option('display.max_columns',30); print(df.round(1).to_string())