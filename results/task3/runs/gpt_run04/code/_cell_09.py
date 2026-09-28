import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('pbmc_processed.h5ad');k='leiden_1.0'
mk=[x for x in ['CD3D','CD3E','IL7R','CCR7','MALAT1','LTB','LST1','LYZ','S100A8','S100A9','FCN1','FCGR3A','MS4A1','CD79A','CD74','TCL1A','NKG7','GNLY','GZMB','GZMK','CCL5','FCER1A','CST3','PPBP','PF4','MKI67','TYMS','TOP2A','HLA-DRA'] if x in ad.var_names]
rows=[]
for c in ad.obs[k].cat.categories:
 x=ad[ad.obs[k]==c]
 rows.append(np.asarray(x[:,mk].X.mean(0)).ravel())
pd.set_option('display.max_columns',40);print(pd.DataFrame(rows,index=ad.obs[k].cat.categories,columns=mk).round(2).to_string())
# check QC means and edges perhaps cluster 2 odd T
print(ad.obs.groupby(k,observed=True)[['total_counts','n_genes_by_counts','pct_counts_mt']].mean().round(1))