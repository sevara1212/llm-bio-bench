import scanpy as sc, numpy as np,pandas as pd
ad=sc.read_h5ad('processed.h5ad');k='leiden_0.4';r=ad.uns['rank_genes_groups']
for c in ['0','4','11','16','17']:
 print('\n',c)
 print([x for x in r['names'][c][:200] if not x.startswith(('RPS','RPL','MT-','HIST'))][:70])
# which clusters each known marker has strongest expression
markers=['CD34','GATA2','AVP','MPO','PRSS57','GYPA','ALAS2','SOX4','DNTT','VPREB1','IGLL1','TCL1A','CCR7','SELL','LEF1','TCF7','MAL','LTB','CCL5','GZMK','TRDC','TRGC1','TRGC2','FCER1A','CD1C','CLEC10A','CD14','S100A8','LILRB1','MS4A7','LILRA4']
for g in markers:
 if g in ad.raw.var_names:
  x=ad.raw[:,g].X.toarray().ravel(); v=pd.Series(x).groupby(ad.obs[k].values).mean()
  print(g, v.idxmax(),round(v.max(),2))