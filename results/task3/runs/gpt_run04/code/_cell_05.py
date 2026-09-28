import scanpy as sc, numpy as np, pandas as pd
sc.settings.verbosity=2
ad=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols',cache=False)
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True,percent_top=None)
ad=ad[(ad.obs.n_genes_by_counts>=200)&(ad.obs.pct_counts_mt<15)].copy()
sc.pp.filter_genes(ad,min_cells=3)
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor='seurat')
print(ad, 'hvg',ad.var.highly_variable.sum())
sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=40)
for r in [.3,.5,.7,1.0]:
 sc.tl.leiden(ad,resolution=r,key_added=f'leiden_{r}',random_state=0)
 print(r,ad.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
# resultant select .5 probably
sc.tl.umap(ad,random_state=0)
ad.write('pbmc_processed.h5ad')