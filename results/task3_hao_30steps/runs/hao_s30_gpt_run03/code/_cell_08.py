import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed.h5ad'); k='leiden_0.6'
genes=['TRAC','CD3D','CD4','CD8A','CCR7','TCF7','MAL','LEF1','IL7R','NKG7','CD247','KLRD1','MS4A1','CD79A','CD14','LYZ','LST1','FCER1A','GZMB','MKI67','TYMS','CD74','GATA1','GATA2','GYPA','GYPC','HBA1','HBA2','HBB','ALAS2','AHSP','PRSS57','KLF1','MPO','KIT','TPSAB1','MS4A2','HDC','FCGR3A']
# raw log data
x=sc.get.obs_df(ad,keys=[k]+[g for g in genes if g in ad.var_names]).groupby(k,observed=True).mean().T
print(x.loc[[g for g in genes if g in x.index]].round(2).to_string())
# percentages positive raw
for c in ['4','13','18']:
 a=ad[ad.obs[k]==c]
 print('\ncluster',c, a.n_obs)
 for g in genes:
  if g in a.var_names:
   z=np.asarray(a.raw[:,g].X).ravel();
   if z.mean()>.2: print(g,round(z.mean(),2),round((z>0).mean(),2),end='; ')
 print()