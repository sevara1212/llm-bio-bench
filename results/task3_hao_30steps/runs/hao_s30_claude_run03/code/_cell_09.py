import scanpy as sc
import pandas as pd

ad = sc.read_h5ad('clustered.h5ad')

sc.tl.rank_genes_groups(ad, 'leiden_05', method='wilcoxon', use_raw=True)

result = ad.uns['rank_genes_groups']
groups = result['names'].dtype.names

top_markers = {}
for g in groups:
    top_markers[g] = list(result['names'][g][:15])

for g in groups:
    print(g, top_markers[g])

ad.write('clustered.h5ad')