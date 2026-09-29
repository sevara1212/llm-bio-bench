import pickle,scanpy as sc,pandas as pd
with open('analysis_results.pkl','rb') as f:a=pickle.load(f)
k='leiden_0.6';sc.tl.rank_genes_groups(a,k,method='wilcoxon',use_raw=False,n_genes=100)
r=sc.get.rank_genes_groups_df(a,group=None)
exclude=('RPL','RPS','MT-','MALAT1','B2M','TMSB10','TMSB4X','ACTB','GAPDH','EEF1A1','FTH1','FTL','IL32','HLA-')
for g in a.obs[k].cat.categories:
 x=r[(r.group==g)&~r.names.str.startswith(exclude)].head(15)
 print(g,':',' '.join(x.names.tolist()))