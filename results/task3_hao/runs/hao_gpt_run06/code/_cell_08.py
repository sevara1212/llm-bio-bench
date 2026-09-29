import scanpy as sc, numpy as np
from scipy.sparse import diags
adata=sc.read_h5ad('pbmc_qc_clustered.h5ad'); X=adata.layers['counts']; s=np.asarray(X.sum(1)).ravel(); Y=diags(1e4/s)@X; Y.data=np.log1p(Y.data); key='leiden_0.5'
genes=['CD4','CD8A','CD8B','TRDC','TRGC1','TRGC2','FOXP3','IL32','NKG7','KLRB1','ZBTB16','SLC4A10','FCGR3A','MS4A7','FCER1A','CD1C','CLEC10A','GZMB','IRF7','IL3RA','LILRA4','MKI67','TOP2A','TCL1A','CD27','CCR7']
ix=[adata.var_names.get_loc(g) for g in genes]
for cl in ['0','1','3','4','6','8','9','10','13','15','16','17','19','20']:
 m=np.asarray(Y[(adata.obs[key]==cl).values][:,ix].mean(0)).ravel()
 print(cl,', '.join(f'{g}={v:.1f}' for g,v in zip(genes,m) if v>.7))