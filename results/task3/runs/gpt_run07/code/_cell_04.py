import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols'); ad.var_names_make_unique()
ad.var['mt']=ad.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False)
print(ad.obs[['n_genes_by_counts','total_counts','pct_counts_mt']].describe(percentiles=[.01,.05,.1,.5,.9,.95,.99]))
print('mt genes',ad.var.mt.sum(), ad.var_names[ad.var.mt][:10].tolist())
for a,b in [(200,500),(200,1000),(200,2500),(200,4000),(500,4000)]: print(a,b,((ad.obs.n_genes_by_counts>=a)&(ad.obs.n_genes_by_counts<=b)&(ad.obs.pct_counts_mt<10)).sum())