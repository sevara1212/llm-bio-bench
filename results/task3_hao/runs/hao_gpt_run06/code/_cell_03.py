import scanpy as sc, numpy as np
adata=sc.read_h5ad('raw_counts.h5ad')
adata.var['mt']=adata.var_names.str.upper().str.startswith('MT-')
adata.var['ribo']=adata.var_names.str.upper().str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt','ribo'],inplace=True,log1p=False)
for c in ['total_counts','n_genes_by_counts','pct_counts_mt','pct_counts_ribo']:
 print('\n',c, np.percentile(adata.obs[c],[0,1,5,10,25,50,75,90,95,99,100]))
print('mt genes',adata.var.mt.sum())
print('cells genes<200',sum(adata.obs.n_genes_by_counts<200),'mt>15',sum(adata.obs.pct_counts_mt>15),'counts >30k',sum(adata.obs.total_counts>30000))