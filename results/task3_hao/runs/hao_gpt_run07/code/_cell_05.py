import scanpy as sc, numpy as np
sc.settings.verbosity=2
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True)
# Standard conservative QC: low complexity, extreme complexity (probable multiplets), high mito
ad=ad[(ad.obs.n_genes_by_counts>=500)&(ad.obs.n_genes_by_counts<=5500)&(ad.obs.pct_counts_mt<=15)].copy()
sc.pp.filter_genes(ad,min_cells=3)
print('post QC',ad)
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4)
sc.pp.log1p(ad)
ad.raw=ad
sc.pp.highly_variable_genes(ad,n_top_genes=3000,flavor='seurat')
print('HVG',ad.var.highly_variable.sum())
sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=40)
for res in [.4,.6,.8,1.0]: sc.tl.leiden(ad,resolution=res,key_added=f'leiden_{res}',flavor='igraph',n_iterations=2)
sc.tl.umap(ad,random_state=0)
for res in [.4,.6,.8,1.0]: print(res,ad.obs[f'leiden_{res}'].value_counts().sort_index().to_dict())
ad.write_h5ad('processed_pbmc.h5ad')