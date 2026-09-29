import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('qc_processed.h5ad'); g='leiden_0.5'
for c in ['4','10','17','0','2','5','15']:
 x=ad[ad.obs[g]==c]
 print('\n',c, x.n_obs)
 for genes in [['MKI67','TOP2A','TUBB','CD3D','NKG7','MS4A1','LYZ','HBB','GATA1','PRSS57','SOX4','DNTT','CD34','TRBC1','TRDC','FCGR3A','GNLY','KLRB1','IFITM3','ISG15','IFIT1','IFIT3']]:
  for z in genes:
   if z in ad.var_names:
    a=x.raw[:,z].X.toarray().ravel(); print(z,round(a.mean(),2),round((a>0).mean(),2),end='; ')
 print()