import pickle,scanpy as sc,pandas as pd
with open('analysis_results.pkl','rb') as f: adata=pickle.load(f)
key='leiden_0.6'
sc.tl.rank_genes_groups(adata,key,method='wilcoxon',use_raw=False,n_genes=30,pts=True)
r=sc.get.rank_genes_groups_df(adata,group=None)
for g in adata.obs[key].cat.categories:
 d=r[r.group==g].head(12)
 print('\n',g,adata.obs[key].value_counts()[g],':',', '.join(d.names.tolist()))
# marker score summary clusters
markers={'CD4_T':['IL7R','LTB','CCR7','MALAT1','LST1'], 'cytotoxic':['NKG7','GNLY','PRF1','GZMB','GZMH'], 'B':['MS4A1','CD79A','CD37','CD74','HLA-DRA'], 'mono':['LYZ','S100A8','S100A9','CTSS','FCN1'], 'DC':['FCER1A','CLEC10A','CD1C','CLEC9A'], 'platelet':['PPBP','PF4','NRGN'], 'IFN':['ISG15','IFIT1','IFIT3','MX1'], 'plasma':['MZB1','JCHAIN','IGHG1','XBP1']}
for n,gs in markers.items():
 valid=[x for x in gs if x in adata.var_names]
 sc.tl.score_genes(adata,valid,score_name=n+'_s',use_raw=False)
print(adata.obs.groupby(key,observed=True)[[x+'_s' for x in markers]].mean().round(2).to_string())