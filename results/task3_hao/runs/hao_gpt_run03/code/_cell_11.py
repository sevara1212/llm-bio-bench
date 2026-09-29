import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad'); k='leiden_0.5'
# raw log expressions marker cluster means
markers=['CD3D','TRAC','NKG7','GNLY','MS4A1','CD79A','LYZ','FCER1A','IL3RA','GZMB','MKI67','TOP2A','CCR7','IL7R','FOXP3','IL2RA','CTLA4','CD14','FCGR3A','CLEC9A','CD34','GATA2','AVP','GYPA','HBB']
df=pd.DataFrame(index=ad.obs[k].cat.categories)
for gene in markers:
 if gene in ad.raw.var_names:
  x=ad.raw[:,gene].X
  import numpy as np
  x=x.toarray().ravel() if hasattr(x,'toarray') else x.ravel()
  df[gene]=pd.Series(x,index=ad.obs[k]).groupby(level=0,observed=True).mean().reindex(df.index)
print(df.round(2).to_string())
# detail cluster4 higher leiden 0.7 correspondence
print('\nCluster4 subclusters'); print(pd.crosstab(ad.obs[k],ad.obs['leiden_0.7']).loc['4'])