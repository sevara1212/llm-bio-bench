import scanpy as sc, numpy as np
adata=sc.read_h5ad('raw_counts.h5ad')
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True,log1p=False)
# conservative filters: sufficient complexity and exclude high mitochondrial / probable doublets
keep=(adata.obs.n_genes_by_counts>=500)&(adata.obs.pct_counts_mt<12)&(adata.obs.n_genes_by_counts<5500)
print('kept',keep.sum(),'removed',(~keep).sum())
print('reasons low/highmt/highgenes', (adata.obs.n_genes_by_counts<500).sum(),(adata.obs.pct_counts_mt>=12).sum(),(adata.obs.n_genes_by_counts>=5500).sum())
adata=adata[keep].copy()
sc.pp.filter_genes(adata,min_cells=3)
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
adata.raw=adata
sc.pp.highly_variable_genes(adata,n_top_genes=3000,flavor='seurat')
print(adata, 'HVG',adata.var.highly_variable.sum())
sc.pp.scale(adata,max_value=10)
sc.tl.pca(adata,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(adata,n_neighbors=15,n_pcs=40)
for r in [.3,.5,.7,.9,1.1]:
 sc.tl.leiden(adata,resolution=r,key_added=f'leiden_{r}',random_state=0)
 print(r,adata.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
adata.write_h5ad('qc_processed.h5ad')