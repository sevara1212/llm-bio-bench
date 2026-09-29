import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed_pbmc.h5ad')
key='leiden_0.6'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',n_genes=15,pts=True)
r=ad.uns['rank_genes_groups']
for cl in r['names'].dtype.names:
 print('\n',cl,ad.obs[key].value_counts()[cl],':', ', '.join(r['names'][cl][:15]))
# marker average table
markers=['IL7R','LTB','CCR7','LST1','S100A8','S100A9','LYZ','FCGR3A','MS4A7','NKG7','GNLY','PRF1','GZMB','CD3D','CD3E','TRAC','MS4A1','CD79A','CD37','CD74','HLA-DRA','GZMK','CCL5','JCHAIN','MZB1','SDC1','PPBP','PF4','FCER1A','CLEC10A','GATA3','KLRD1','TYROBP','IFIT1','ISG15','MKI67']
print('\nMEANS')
print(pd.DataFrame(ad[:,markers].X.toarray(),index=ad.obs[key]).groupby(level=0,observed=True).mean().round(1).to_string())