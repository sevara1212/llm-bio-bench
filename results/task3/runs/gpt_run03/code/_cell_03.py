import scanpy as sc, numpy as np, pandas as pd
p='filtered_gene_bc_matrices/hg19'
adata=sc.read_10x_mtx(p, var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
print(adata)
print(adata.var.head())
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True,log1p=False)
print(adata.obs[['n_genes_by_counts','total_counts','pct_counts_mt']].describe(percentiles=[.01,.05,.1,.25,.5,.75,.9,.95,.99]))
print('genes',adata.var_names[:20].tolist(), 'mt genes',adata.var.mt.sum())
# distributions / quantile upper
for col in ['n_genes_by_counts','total_counts','pct_counts_mt']:
 print(col, np.quantile(adata.obs[col],[0,.001,.01,.05,.95,.99,.999,1]))