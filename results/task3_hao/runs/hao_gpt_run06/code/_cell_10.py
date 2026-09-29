import scanpy as sc, numpy as np
raw=sc.read_h5ad('raw_counts.h5ad'); raw.var['mt']=raw.var_names.str.startswith('MT-');sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,log1p=False)
bad=raw[(raw.obs.n_genes_by_counts<500)|(raw.obs.pct_counts_mt>=15)]
print(bad.obs[['total_counts','n_genes_by_counts','pct_counts_mt']])
for i in range(bad.n_obs):
 x=bad.X[i].toarray().ravel() if hasattr(bad.X,'toarray') else bad.X[i]
 inds=np.argsort(x)[-20:][::-1]
 print('\n',bad.obs_names[i]);print([(bad.var_names[j],x[j]) for j in inds[:15]])