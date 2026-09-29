import scanpy as sc, numpy as np, pandas as pd, gc
ad=sc.read_h5ad('raw_counts.h5ad'); ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-');sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True)
ad=ad[(ad.obs.n_genes_by_counts>=500)&(ad.obs.n_genes_by_counts<=5500)&(ad.obs.pct_counts_mt<12)].copy();sc.pp.filter_genes(ad,min_cells=3)
sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad); logX=ad.X.copy()
sc.pp.highly_variable_genes(ad,n_top_genes=3000,flavor='seurat'); sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50,use_highly_variable=True,svd_solver='arpack');sc.pp.neighbors(ad,n_neighbors=15,n_pcs=30);sc.tl.umap(ad,random_state=0)
for r in [.4,.6,.8,1.0]:
 sc.tl.leiden(ad,resolution=r,key_added=f'leiden_{r}',random_state=0)
 print('resolution',r, ad.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
ad.X=logX
# .6 chooses expected pbmc subtypes
sc.tl.rank_genes_groups(ad,groupby='leiden_0.6',method='wilcoxon',use_raw=False,n_genes=30)
for g in ad.obs['leiden_0.6'].cat.categories:
 print('\nCLUSTER',g,'n=',sum(ad.obs['leiden_0.6']==g))
 print(', '.join(ad.uns['rank_genes_groups']['names'][g][:15]))
# lightweight stored artifacts marker table, metadata and barcode clusters
pd.DataFrame(ad.uns['rank_genes_groups']['names']).to_csv('cluster_markers_top30.csv',index=False)
ad.obs[['leiden_0.6']].to_csv('cluster_assignments_qc.csv')