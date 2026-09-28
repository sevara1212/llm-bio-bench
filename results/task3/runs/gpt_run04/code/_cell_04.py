import scanpy as sc, numpy as np, pandas as pd
sc.settings.verbosity=2
adata=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols',cache=False)
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
adata.var['ribo']=adata.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt','ribo'],inplace=True,percent_top=None,log1p=False)
print(adata.obs[['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']].describe(percentiles=[.01,.05,.5,.95,.99]))
print('mt genes',adata.var.mt.sum())
# Criteria conservative: cells at least 200 genes, mt <=15; genes 3 cell.
keep=(adata.obs.n_genes_by_counts>=200)&(adata.obs.pct_counts_mt<15)
print('keep',keep.sum(),'remove',(~keep).sum())
print(adata.obs.loc[~keep,['total_counts','n_genes_by_counts','pct_counts_mt']].sort_values('pct_counts_mt').tail())