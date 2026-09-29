import scanpy as sc, pandas as pd
raw=sc.read_h5ad('raw_counts.h5ad'); raw.var['mt']=raw.var_names.str.upper().str.startswith('MT-');sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,log1p=False)
ok=(raw.obs.n_genes_by_counts>=500)&(raw.obs.pct_counts_mt<15)
print(raw.obs.loc[~ok,['n_genes_by_counts','total_counts','pct_counts_mt']])
# show marker expression for excluded
for bc in raw.obs_names[~ok]:
 x=raw[bc,['CD3D','NKG7','MS4A1','LYZ','PPBP','HBB','MKI67']].X.toarray().ravel()
 print(bc,dict(zip(['CD3D','NKG7','MS4A1','LYZ','PPBP','HBB','MKI67'],x)))