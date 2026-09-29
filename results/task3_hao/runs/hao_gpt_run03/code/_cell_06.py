import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('raw_counts.h5ad')
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True,log1p=False)
# retain robust cells; all pass gene cutoff almost
ad=ad[(ad.obs.n_genes_by_counts>=500)&(ad.obs.pct_counts_mt<15)].copy()
sc.pp.filter_genes(ad,min_cells=3)
ad.obs['batch']=ad.obs_names.str.split('_').str[0].astype('category')
# Preserve raw counts then normalize
ad.layers['counts']=ad.X.copy()
sc.pp.normalize_total(ad,target_sum=1e4)
sc.pp.log1p(ad)
ad.raw=ad
# batch-aware HVGs
sc.pp.highly_variable_genes(ad,flavor='seurat',n_top_genes=3000,batch_key='batch')
print(ad, 'HVG',ad.var.highly_variable.sum())
# regression potential makes clusters better
sc.pp.regress_out(ad,['total_counts','pct_counts_mt'])
sc.pp.scale(ad,max_value=10)
sc.tl.pca(ad,n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=40)
sc.tl.umap(ad,random_state=42)
for r in [.3,.5,.7,.9,1.1]: sc.tl.leiden(ad,resolution=r,key_added=f'leiden_{r}',random_state=42)
ad.write_h5ad('processed.h5ad')
for r in [.3,.5,.7,.9,1.1]: print(r,ad.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())