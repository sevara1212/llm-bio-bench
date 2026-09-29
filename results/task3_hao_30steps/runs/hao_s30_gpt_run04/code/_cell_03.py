import scanpy as sc, numpy as np, pandas as pd
adata=sc.read_h5ad('raw_counts.h5ad')
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
adata.var['ribo']=adata.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt','ribo'],inplace=True,log1p=False)
for x in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 print(x, np.percentile(adata.obs[x],[0,1,5,25,50,75,95,99,100]))
print('MT genes',adata.var.mt.sum())
print(adata.obs.sort_values('pct_counts_mt',ascending=False)[['total_counts','n_genes_by_counts','pct_counts_mt']].head(10))
# genes seen & cell counts
print('cells/gene percentiles',np.percentile(adata.var.n_cells_by_counts,[0,1,25,50,75,95,99,100]))