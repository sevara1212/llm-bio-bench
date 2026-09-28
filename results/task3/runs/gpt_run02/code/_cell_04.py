import scanpy as sc, numpy as np, pandas as pd
p='filtered_gene_bc_matrices/hg19'; adata=sc.read_10x_mtx(p,var_names='gene_symbols')
adata.var_names_make_unique(); adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True)
# retain cellular barcodes; low quality cutoff on mitochondrial fraction, plus expressed genes
adata=adata[(adata.obs.n_genes_by_counts>=200)&(adata.obs.pct_counts_mt<10)].copy()
adata.layers['counts']=adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=2000)
print('post QC',adata.shape,'HVG',adata.var.highly_variable.sum())
sc.pp.scale(adata,max_value=10); sc.tl.pca(adata,n_comps=50,use_highly_variable=True,svd_solver='arpack')
sc.pp.neighbors(adata,n_neighbors=15,n_pcs=30)
sc.tl.umap(adata,random_state=0)
for r in [.4,.5,.6,.8,1.0]:
 sc.tl.leiden(adata,resolution=r,key_added=f'leiden_{r}',random_state=0)
 print(r,adata.obs[f'leiden_{r}'].value_counts().sort_index().to_dict())
# save preliminary
adata.write('pbmc_processed.h5ad')