import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_counts.h5ad'); ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True); ad=ad[ad.obs.pct_counts_mt<15].copy(); sc.pp.filter_genes(ad,min_cells=3)
sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad);ad.raw=ad
sc.pp.highly_variable_genes(ad,n_top_genes=3000,flavor='seurat')
# scale only hvg copy, run PCA and transfer pca
h=ad[:,ad.var.highly_variable].copy();sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=50,svd_solver='arpack')
ad.obsm['X_pca']=h.obsm['X_pca']; ad.varm['PCs']=np.zeros((ad.n_vars,50));ad.varm['PCs'][ad.var.highly_variable.values]=h.varm['PCs']
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=40);sc.tl.umap(ad,random_state=42)
for r in [.3,.4,.5,.6,.7,.8,1]:
 sc.tl.leiden(ad,resolution=r,key_added=f'leiden_{r}',random_state=42)
 print(r,ad.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
# use .5
ad.obs['cluster']=ad.obs['leiden_0.5'].astype(str)
sc.tl.rank_genes_groups(ad,'cluster',method='wilcoxon',use_raw=True,n_genes=20)
for cl in sorted(ad.obs.cluster.unique(),key=int):
 d=sc.get.rank_genes_groups_df(ad,group=cl).head(12)
 print('\nCLUSTER',cl,'N',sum(ad.obs.cluster==cl));print(', '.join(d.names.tolist()))
# save lean, no dense h data
import pickle
with open('clustered.pkl','wb') as f:pickle.dump(ad,f)