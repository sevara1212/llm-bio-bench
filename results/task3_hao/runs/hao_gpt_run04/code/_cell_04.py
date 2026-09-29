import scanpy as sc, numpy as np
adata=sc.read_h5ad('raw_counts.h5ad')
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True, log1p=False)
# QC: low-complexity / high-mitochondrial cells; retain high RNA cells since PBMC activation/types vary
keep=(adata.obs.n_genes_by_counts>=500)&(adata.obs.pct_counts_mt<12)
print('kept',keep.sum(),'of',adata.n_obs,'removed',(~keep).sum())
print('removed reasons', (adata.obs.n_genes_by_counts<500).sum(), (adata.obs.pct_counts_mt>=12).sum())
adata=adata[keep].copy()
adata.obs['batch']=adata.obs_names.str.split('_').str[0].astype('category')
adata.layers['counts']=adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
try:
 sc.pp.highly_variable_genes(adata,n_top_genes=3000,flavor='seurat_v3',layer='counts',batch_key='batch')
except Exception as e:
 print('v3 fail',repr(e)); sc.pp.highly_variable_genes(adata,n_top_genes=3000,flavor='seurat')
print('HVG',adata.var.highly_variable.sum())
sc.pp.scale(adata,max_value=10)
sc.tl.pca(adata,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(adata,n_neighbors=15,n_pcs=40)
sc.tl.umap(adata,random_state=0)
for res in [.4,.6,.8,1.0]:
 sc.tl.leiden(adata,resolution=res,key_added=f'leiden_{res}',random_state=0)
 print(res,adata.obs[f'leiden_{res}'].value_counts().sort_index().to_dict())
adata.write_h5ad('processed.h5ad')