import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True,log1p=False)
# retain robust cells; remove 3 cells with extreme mitochondrial percentage
ad=ad[(ad.obs.n_genes_by_counts>=500)&(ad.obs.pct_counts_mt<15)].copy()
print('retained',ad.shape, 'removed',7841-ad.n_obs)
# Preserve counts in raw then normalize
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4)
sc.pp.log1p(ad)
ad.raw=ad
sc.pp.highly_variable_genes(ad,n_top_genes=3000,flavor='seurat')
print('HVG',ad.var.highly_variable.sum())
sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50,use_highly_variable=True,svd_solver='arpack',random_state=0)
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=40,random_state=0)
sc.tl.umap(ad,random_state=0)
for res in [.4,.5,.6,.7,.8,1.0]:
 sc.tl.leiden(ad,resolution=res,key_added=f'leiden_{res}',random_state=0)
 print(res,ad.obs[f'leiden_{res}'].value_counts().sort_index().to_dict())
ad.write_h5ad('processed.h5ad')