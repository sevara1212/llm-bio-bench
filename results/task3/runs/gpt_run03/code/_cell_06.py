import scanpy as sc, numpy as np, pandas as pd
A=sc.read_h5ad('pbmc_processed_r09.h5ad'); cl=A.obs.leiden
# reconstruct normalized expression from counts layer
B=A.copy();B.X=B.layers['counts'].copy();sc.pp.normalize_total(B,target_sum=1e4);sc.pp.log1p(B)
genes=['CD3D','CD3E','IL7R','LTB','CCR7','MALAT1','TRBC2','LST1','LYZ','S100A8','S100A9','FCGR3A','MS4A7','NKG7','GNLY','PRF1','GZMB','GZMK','CCL5','CD79A','MS4A1','CD74','HLA-DRA','TCL1A','FCER1A','CLEC10A','PPBP','PF4','FCGR3B','JCHAIN','MZB1']
inds=[B.var_names.get_loc(g) for g in genes if g in B.var_names]
df=pd.DataFrame(B.X[:,inds].toarray(),index=cl,columns=[B.var_names[i] for i in inds]).groupby(level=0,observed=True).mean().T
print(df.round(2).to_string())
# rank based B
sc.tl.rank_genes_groups(B,'leiden',method='wilcoxon',n_genes=25)
for c in sorted(cl.unique(),key=int): print(c, ', '.join(sc.get.rank_genes_groups_df(B,group=c).head(15).names))