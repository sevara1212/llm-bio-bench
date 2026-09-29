import scanpy as sc,pandas as pd,numpy as np
ad=sc.read_h5ad('raw_counts.h5ad'); cl=pd.read_csv('cluster_assignments_qc.csv',index_col=0).iloc[:,0]; ad=ad[cl.index].copy();sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad)
genes=['CD3D','TRAC','IL7R','CCR7','LTB','MALAT1','NKG7','GNLY','KLRD1','MS4A1','CD79A','CD74','CD14','LYZ','S100A8','FCGR3A','LST1','FCER1A','CLEC10A','CLEC9A','GZMB','GZMK','MKI67','TYMS','HBA1','HBB','PRSS57','PF4','PPBP','JCHAIN','MZB1','TRDC','KLRB1','IFITM3']
X=ad[:,genes].X
for c in cl.astype(str).sort_values().unique():
 ix=np.where(cl.astype(str).values==c)[0]
 y=np.asarray(X[ix].mean(axis=0)).ravel()
 print(c, ' | '.join(f'{g}:{v:.1f}' for g,v in zip(genes,y) if v>.5))