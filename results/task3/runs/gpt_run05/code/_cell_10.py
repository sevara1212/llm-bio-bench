import scanpy as sc, numpy as np
p='filtered_gene_bc_matrices/hg19/'
a=sc.read_10x_mtx(p,var_names='gene_symbols');a.var_names_make_unique();a.var['mt']=a.var_names.str.startswith('MT-');sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
x=a.obs[(a.obs.n_genes_by_counts<200)|(a.obs.pct_counts_mt>=15)]
print(x[['n_genes_by_counts','total_counts','pct_counts_mt']])
for b in x.index:
 i=a.obs_names.get_loc(b); inds=a.X[i].indices; dat=a.X[i].data; order=np.argsort(dat)[::-1][:20]
 print('\n',b, list(zip(a.var_names[inds[order]],dat[order])))