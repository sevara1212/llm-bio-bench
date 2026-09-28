import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('pbmc_qc_clustered.h5ad')
# restore log-normalized expression for DE
ad.X=ad.layers['counts'].copy(); sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
sc.tl.rank_genes_groups(ad,'leiden',method='wilcoxon',n_genes=30)
r=ad.uns['rank_genes_groups']
for c in ad.obs.leiden.cat.categories:
 print('\nCLUSTER',c,'n=',sum(ad.obs.leiden==c))
 print(', '.join(r['names'][c][:20]))
# inspect canonical markers cluster mean log expr
markers=['IL7R','LTB','CCR7','LST1','S100A8','S100A9','LYZ','FCGR3A','MS4A7','NKG7','GNLY','PRF1','GZMB','CD3D','TRAC','CD79A','MS4A1','CD37','CD74','HLA-DRA','GZMK','MALAT1','PPBP','PF4','FCER1A','CST3']
D=pd.DataFrame(index=ad.obs.leiden.cat.categories)
for x in markers:
 D[x]=np.asarray(ad[:,x].X.mean(0)).ravel()[0] if x in ad.var_names else np.nan
print('\nMEAN lognorm'); print(D.round(2).to_string())
ad.write('pbmc_qc_clustered.h5ad')