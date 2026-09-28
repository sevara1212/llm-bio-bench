import scanpy as sc,numpy as np,pandas as pd
ad=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols'); ad.var_names_make_unique()
ad.var['mt']=ad.var_names.str.startswith('MT-');ad.var['ribo']=ad.var_names.str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(ad,qc_vars=['mt','ribo'],inplace=True,log1p=False)
for x in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 q=np.percentile(ad.obs[x],[0,1,5,25,50,75,95,99,100]);print(x, np.round(q,2).tolist())
print('threshold failures', (ad.obs.n_genes_by_counts<200).sum(),(ad.obs.pct_counts_mt>15).sum(),(ad.obs.total_counts>30000).sum())
print(ad.obs.sort_values('pct_counts_mt',ascending=False)[['total_counts','n_genes_by_counts','pct_counts_mt']].head(10))