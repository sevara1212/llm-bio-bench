import scanpy as sc, pandas as pd, numpy as np
ad=sc.read_h5ad('processed_markers.h5ad'); key='leiden_0.5'
genes=['CD3D','CD3E','TRAC','CCR7','LTB','MAL','LEF1','TCF7','IL7R','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A7','NKG7','GNLY','KLRD1','GZMB','MS4A1','CD79A','CD74','HLA-DRA','JCHAIN','MZB1','GZMK','KLRB1','TRDC','TRGC1','MKI67','TOP2A','CD34','GATA2','PRSS57','HBB','PPBP','PF4','CLEC9A','FCER1A','LILRA4','GZMB']
# raw is log-normalized unscaled
x=ad.raw[:,genes].X
if hasattr(x,'toarray'):x=x.toarray()
df=pd.DataFrame(x,index=ad.obs[key].astype(str),columns=genes).groupby(level=0,observed=True).mean().T
print(df.round(2).to_string())
# Find cl0 specific excluding RPS RPL
r=ad.uns['rank_genes_groups']['names']
print('\ncluster 0 expanded', [x for x in r['0'] if not x.startswith(('RP','MT-'))][:40])