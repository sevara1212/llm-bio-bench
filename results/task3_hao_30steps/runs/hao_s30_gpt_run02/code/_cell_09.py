import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed_markers.h5ad'); key='leiden_0.6'
genes=['TRAC','CD3D','NKG7','KLRD1','MS4A1','CD79A','LYZ','LST1','MPO','HBB','GYPA','KLRB1','IL7R','CD8A','MKI67','TOP2A','TYMS','PTPRC']
X=ad.raw[:,genes].X
if hasattr(X,'toarray'): X=X.toarray()
df=pd.DataFrame(X,index=ad.obs[key],columns=genes).groupby(level=0,observed=True).mean()
print(df.round(2).to_string())