import scanpy as sc, pandas as pd, numpy as np
adata=sc.read_h5ad('qc_processed.h5ad'); grp='leiden_0.5'
genes=['CD3D','CD3E','TRAC','IL7R','LTB','CCR7','MALAT1','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A7','IFITM3','FCER1A','CLEC10A','CD1C','CLEC9A','GZMB','NKG7','GNLY','KLRD1','KLRB1','MS4A1','CD79A','CD74','HLA-DRA','MZB1','JCHAIN','IGKC','CD27','GZMK','TCF7','MKI67','TYMS','TRDC','FCGR3A','PPBP','PF4','HBB','LILRA4','GZMB']
D=pd.DataFrame(adata.raw[:,genes].X.toarray(),index=adata.obs[grp],columns=genes).groupby(level=0,observed=True).mean()
pd.set_option('display.max_columns',100); print(D.round(1).T.to_string())
# full top non ribosomal
sc.tl.rank_genes_groups(adata,groupby=grp,method='wilcoxon',use_raw=True,n_genes=100)
for c in adata.obs[grp].cat.categories:
 gg=[x for x in adata.uns['rank_genes_groups']['names'][c] if not x.startswith(('RPL','RPS','MT-'))][:12]
 print(c,':',','.join(gg))