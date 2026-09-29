import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_counts.h5ad'); ad.var['mt']=ad.var_names.str.startswith('MT-'); sc.pp.calculate_qc_metrics(ad,qc_vars=['mt'],inplace=True)
for c in ['pct_counts_mt','n_genes_by_counts','total_counts']:
 print(c, np.percentile(ad.obs[c],[0,1,5,50,90,95,99,99.5,99.9,100]))
for t in [10,15,20,25]:print(t,(ad.obs.pct_counts_mt>t).sum())
print('cells hbb>30', ((ad[:,ad.var_names=='HBB'].X.toarray().ravel()/ad.obs.total_counts*100)>30).sum())