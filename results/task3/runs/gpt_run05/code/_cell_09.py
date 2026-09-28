import scanpy as sc, numpy as np, pandas as pd
A=sc.read_h5ad('pbmc_processed_prelim.h5ad'); cl='leiden_1.0'
# fraction detected and mean counts by group for discriminating genes
markers=['CCR7','LTB','MALAT1','IL7R','LST1','CD3D','CD3E','TRAC','CCR10','MKI67','NKG7','GNLY','MS4A1','CD79A','CD14','LYZ','FCGR3A','FCER1A','CLEC10A','PPBP','PF4','S100A8','S100A9','GZMB','CCL5','CD4','CD8A']
X=A.layers['counts']; genes=list(A.var_names)
for g in sorted(A.obs[cl].unique(),key=int):
 ix=np.where(A.obs[cl].values==g)[0]
 rows=[]
 for m in markers:
  j=genes.index(m) if m in genes else None
  if j is not None:
   v=X[ix,j].toarray().ravel() if hasattr(X[ix,j],'toarray') else X[ix,j]
   if (v>0).mean()>.1: rows.append(f'{m}:{(v>0).mean():.2f}/{v.mean():.1f}')
 print(g, ' '.join(rows))