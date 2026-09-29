import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.startswith('MT-'); sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True)
# conservative QC: retain real PBMCs but discard low quality/high mitochondrial
ad=ad[(ad.obs.n_genes_by_counts>=500)&(ad.obs.pct_counts_mt<15)].copy()
print('retained',ad.shape)
# Remove genes absent/rare; retain raw counts in layer
sc.pp.filter_genes(ad,min_cells=3)
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
ad.raw=ad
sc.pp.highly_variable_genes(ad,n_top_genes=3000,flavor='seurat')
print('HVG',ad.var.highly_variable.sum())
sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50,use_highly_variable=True,svd_solver='arpack',random_state=0)
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=35,random_state=0)
for r in [.4,.6,.8,1.0]: sc.tl.leiden(ad,resolution=r,key_added=f'leiden_{r}',random_state=0)
for r in [.4,.6,.8,1.0]: print(r,ad.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
ad.write('processed.h5ad')