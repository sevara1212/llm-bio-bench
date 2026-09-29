import scanpy as sc, pandas as pd, numpy as np
adata=sc.read_h5ad('pbmc_qc_clustered.h5ad'); key='leiden_0.5'
# Use lognormalized X (scaled unfortunately X scaled); layer counts raw. recreate log norm in separate raw counts aggregate
X=adata.layers['counts']
# normalized log cell-by-gene sparse
from scipy.sparse import diags
s=np.asarray(X.sum(1)).ravel(); Y=diags(1e4/s)@X; Y.data=np.log1p(Y.data)
genes=['CD3D','TRAC','IL7R','CCR7','LTB','MAL','LEF1','TCF7','LST1','LYZ','S100A8','MS4A1','CD79A','NKG7','GNLY','GZMB','GZMK','KLRD1','HLA-DRA','FCER1A','CLEC9A','MKI67','STMN1','TYMS','HBB','PRSS57','SOX4','CD74']
ix=[adata.var_names.get_loc(g) for g in genes]
for cl in adata.obs[key].cat.categories:
 mask=(adata.obs[key]==cl).values
 means=np.asarray(Y[mask][:,ix].mean(0)).ravel()
 top=sorted(zip(genes,means),key=lambda z:-z[1])[:8]
 print(cl, ' '.join(f'{g}:{v:.1f}' for g,v in top))