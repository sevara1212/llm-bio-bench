import scanpy as sc, pandas as pd, numpy as np
p='filtered_gene_bc_matrices/hg19/'
a=sc.read_10x_mtx(p,var_names='gene_symbols'); a.var_names_make_unique()
a.var['mt']=a.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
# conservative QC: low genes/debris or extreme mito
keep=(a.obs.n_genes_by_counts>=200)&(a.obs.pct_counts_mt<15)
print('retained',keep.sum(),'of',a.n_obs)
a=a[keep].copy()
a.layers['counts']=a.X.copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat')
print('HVG',a.var.highly_variable.sum())
sc.pp.scale(a,max_value=10)
sc.tl.pca(a,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(a,n_neighbors=15,n_pcs=30)
for r in [.4,.5,.6,.7,.8,1.0]:
 sc.tl.leiden(a,resolution=r,key_added='leiden_'+str(r),random_state=0)
 print(r,a.obs['leiden_'+str(r)].value_counts().sort_index().to_dict())
# choose .6 initially
sc.tl.leiden(a,resolution=.6,key_added='leiden',random_state=0)
sc.tl.rank_genes_groups(a,'leiden',method='wilcoxon',use_raw=False,n_genes=15)
print(sc.get.rank_genes_groups_df(a,group=None).groupby('group').head(12).to_string(index=False))
a.write('pbmc_processed_prelim.h5ad')