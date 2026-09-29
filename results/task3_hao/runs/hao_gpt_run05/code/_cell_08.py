import scanpy as sc,pandas as pd,numpy as np
ad=sc.read_h5ad('processed.h5ad'); key='leiden_0.4'
genes=['CD3D','CD3E','TRAC','IL7R','CCR7','LTB','MALAT1','NKG7','GNLY','MS4A1','CD79A','CD37','CD14','LYZ','S100A8','FCGR3A','LST1','FCER1A','CLEC9A','GZMB','MKI67','TOP2A','CDK1','CD34','GATA2','HLA-DRA','IFITM3','TCL1A','IGHM','CD27','TRDC','FOXP3','IL32']
X=ad.raw[:,genes].X
if hasattr(X,'toarray'): X=X.toarray()
df=pd.DataFrame(X,index=ad.obs[key],columns=genes).groupby(level=0,observed=True).mean()
pd.set_option('display.max_columns',100); print(df.round(2).T)
# non ribo top using score ranking
r=ad.uns['rank_genes_groups']
for c in ['0','4','11']:
 ns=r['names'][c]; print(c,[g for g in ns if not g.startswith(('RPS','RPL','MT-'))][:30])