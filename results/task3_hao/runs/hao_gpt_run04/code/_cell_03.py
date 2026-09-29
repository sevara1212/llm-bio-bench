import scanpy as sc, numpy as np, pandas as pd
adata=sc.read_h5ad('raw_counts.h5ad')
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
adata.var['ribo']=adata.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt','ribo'],inplace=True, log1p=False)
for c in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 print(c, adata.obs[c].quantile([0,.01,.05,.1,.25,.5,.75,.9,.95,.99,1]).round(2).to_dict())
print('zero genes',sum(adata.obs.n_genes_by_counts==0),'mt genes',adata.var.mt.sum())
print(adata.obs[['total_counts','n_genes_by_counts','pct_counts_mt']].describe())
print('prefix counts',adata.obs_names.str.split('_').str[0].value_counts())