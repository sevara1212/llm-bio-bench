import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed_pbmc.h5ad'); key='leiden_0.6'
# calculate mean raw, display primary lineage genes/cell-cycle for cl 0/1 plus all
z=ad.raw.to_adata(); z.obs[key]=ad.obs[key]
genes=['CD3D','TRAC','CD247','CD4','CD8A','CD8B','NKG7','KLRD1','MS4A1','CD79A','CD74','HLA-DRA','LYZ','LST1','S100A8','PF4','PPBP','HBB','HBA1','MPO','ELANE','PRSS57','GATA1','KLF1','SOX4','CD34','MKI67','TOP2A','TYMS','STMN1','TCL1A','IGHM','FCGR3A']
df=pd.DataFrame(z[:,genes].X.toarray(),index=z.obs[key]).groupby(level=0,observed=True).mean()
print(df.loc[['0','1','6','8','11','12','15','16','17','18','20']].round(2).to_string())
# score 0 individual counts of highest lineage expression