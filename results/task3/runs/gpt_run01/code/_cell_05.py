import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols');ad.var_names_make_unique()
ad.var['mt']=ad.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True,log1p=False)
# retain all cells meeting 200 genes, and mt<=15%; genes expressed in >=3
ad=ad[(ad.obs.n_genes_by_counts>=200)&(ad.obs.pct_counts_mt<=15),:].copy(); sc.pp.filter_genes(ad,min_cells=3)
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad,n_top_genes=2000,flavor='seurat')
sc.tl.pca(ad,use_highly_variable=True,n_comps=50,svd_solver='arpack')
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=40)
sc.tl.umap(ad,random_state=0)
for r in [0.4,0.5,0.6,0.7,0.8,1.0]:
 sc.tl.leiden(ad,resolution=r,key_added=f'leiden_{r}',random_state=0)
 print(r,ad.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
# choose .6
ad.obs['cluster']=ad.obs['leiden_0.6'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=False,n_genes=30)
for cl in sorted(ad.obs.cluster.unique(),key=int):
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(15)
 print('\nCLUSTER',cl,'n=',(ad.obs.cluster==cl).sum());print(', '.join(d.names.tolist()))
ad.write('pbmc_processed.h5ad')
print(ad)