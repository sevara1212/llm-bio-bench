import scanpy as sc, pandas as pd, numpy as np
p='filtered_gene_bc_matrices/hg19/'
adata=sc.read_10x_mtx(p,var_names='gene_symbols', cache=False)
adata.var_names_make_unique()
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
adata.var['ribo']=adata.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False)
print(adata)
print(adata.obs[['total_counts','n_genes_by_counts','pct_counts_mt']].describe(percentiles=[.01,.05,.5,.95,.99]))
print('mt genes',adata.var.mt.sum())
# number by filtering
for g, m, lo, hi in [(200,20,0,1e9),(300,15,500,15000),(400,15,500,15000)]:
 k=(adata.obs.n_genes_by_counts>=g)&(adata.obs.pct_counts_mt<m)&(adata.obs.total_counts>=lo)&(adata.obs.total_counts<=hi)
 print(g,m,lo,hi,k.sum(), 'removed',(~k).sum())