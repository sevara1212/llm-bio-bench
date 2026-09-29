import scanpy as sc, pandas as pd
ad=sc.read_h5ad('processed.h5ad')
key='leiden_0.5'
sc.tl.rank_genes_groups(ad,key,method='wilcoxon',use_raw=True,n_genes=15,pts=True)
r=ad.uns['rank_genes_groups']
groups=r['names'].dtype.names
for g in groups:
 genes=list(r['names'][g][:12]); scores=list(r['scores'][g][:3]); n=(ad.obs[key].astype(str)==g).sum()
 print(g,n,', '.join(genes))
ad.write_h5ad('processed_markers.h5ad')