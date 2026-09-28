import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols');ad.var_names_make_unique();ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False)
# retain cells passing conventional low-quality / high-complexity / mito criteria
ad=ad[(ad.obs.n_genes_by_counts>=200)&(ad.obs.n_genes_by_counts<=2500)&(ad.obs.pct_counts_mt<5)].copy()
print('retained',ad.shape, 'removed',2700-ad.n_obs)
sc.pp.filter_genes(ad,min_cells=3)
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad)
ad.raw=ad
sc.pp.highly_variable_genes(ad,flavor='seurat',n_top_genes=2000)
print('genes after expression filter',ad.n_vars,'HVG',ad.var.highly_variable.sum())
ad=ad[:,ad.var.highly_variable].copy()
sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50,svd_solver='arpack'); sc.pp.neighbors(ad,n_neighbors=10,n_pcs=30)
for r in [.3,.4,.5,.6,.7,.8,1.0]:
 sc.tl.leiden(ad,resolution=r,key_added='leiden_'+str(r),random_state=0)
 print(r,ad.obs['leiden_'+str(r)].value_counts().sort_index().to_dict())
# save
ad.write('pbmc_processed_hvg.h5ad')