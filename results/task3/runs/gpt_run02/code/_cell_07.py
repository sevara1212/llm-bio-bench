import scanpy as sc, numpy as np, pandas as pd
ad=sc.read_h5ad('pbmc_processed.h5ad')
print(ad.layers['counts'].min(),ad.layers['counts'].max(),ad.X.min(),ad.X.max())
# excluded QC cells and canonical raw counts
raw=sc.read_10x_mtx('filtered_gene_bc_matrices/hg19',var_names='gene_symbols');raw.var_names_make_unique();raw.var['mt']=raw.var_names.str.startswith('MT-');sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True)
bad=raw[(raw.obs.n_genes_by_counts<200)|(raw.obs.pct_counts_mt>=10)].copy();print(bad.obs[['n_genes_by_counts','total_counts','pct_counts_mt']])
for b in bad.obs_names:
 x=raw[b,:]
 ix=np.asarray(x.X.toarray()).ravel().argsort()[-20:][::-1]
 print('\n',b, list(zip(x.var_names[ix],np.asarray(x.X.toarray()).ravel()[ix])))