import os,glob
for f in glob.glob('processed.h5ad*'):
 print('remove',f,os.path.getsize(f));os.remove(f)
print('free');os.system('df -h .')
import scanpy as sc
adata=sc.read_h5ad('raw_counts.h5ad')
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True,log1p=False)
adata=adata[(adata.obs.n_genes_by_counts>=500)&(adata.obs.pct_counts_mt<12)].copy()
adata.obs['batch']=adata.obs_names.str.split('_').str[0].astype('category')
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=3000,flavor='seurat')
print('kept / HVGs',adata.n_obs,adata.var.highly_variable.sum())
sc.tl.pca(adata,n_comps=50,mask_var='highly_variable',zero_center=False,svd_solver='arpack')
sc.pp.neighbors(adata,n_neighbors=15,n_pcs=40)
sc.tl.umap(adata,random_state=0)
for res in [.4,.6,.8,1.0]:
 sc.tl.leiden(adata,resolution=res,key_added=f'leiden_{res}',random_state=0)
 print('RES',res,adata.obs[f'leiden_{res}'].value_counts().sort_index().to_dict())
# save lean to pickle? output only no need
import pickle
with open('analysis_results.pkl','wb') as f: pickle.dump(adata,f,protocol=4)
print('pickle MB',os.path.getsize('analysis_results.pkl')/1e6)