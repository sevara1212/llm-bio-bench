import scanpy as sc, numpy as np, pandas as pd
p='filtered_gene_bc_matrices/hg19'; a=sc.read_10x_mtx(p,var_names='gene_symbols'); a.var_names_make_unique()
a.var['mt']=a.var_names.str.upper().str.startswith('MT-'); sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],inplace=True)
# conservative cell QC: retain viable cells; high upper gene cutoff excludes likely multiplets
a=a[(a.obs.n_genes_by_counts>=300)&(a.obs.n_genes_by_counts<=2500)&(a.obs.pct_counts_mt<10)].copy()
print('retained',a.n_obs,'/',2700)
a.layers['counts']=a.X.copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat')
print('HVG',a.var.highly_variable.sum())
sc.pp.scale(a,max_value=10); sc.tl.pca(a,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(a,n_neighbors=15,n_pcs=30); sc.tl.leiden(a,resolution=.5,key_added='leiden',random_state=0)
sc.tl.rank_genes_groups(a,'leiden',method='wilcoxon',use_raw=False,n_genes=30)
print(a.obs.leiden.value_counts().sort_index())
for cl in sorted(a.obs.leiden.unique(),key=int):
 d=sc.get.rank_genes_groups_df(a,group=cl).head(12)
 print('\nCLUSTER',cl,'n=',sum(a.obs.leiden==cl)); print(', '.join(d.names.tolist()))
# save processed for inspection
# Note raw absent but layers counts.
a.write('pbmc_processed.h5ad')