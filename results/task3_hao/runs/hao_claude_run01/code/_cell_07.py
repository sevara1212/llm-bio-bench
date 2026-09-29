import scanpy as sc
ad = sc.read_h5ad('processed.h5ad')

sc.tl.rank_genes_groups(ad, 'leiden', method='wilcoxon')

import pandas as pd
result = ad.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = {}
for g in groups:
    top_genes[g] = list(result['names'][g][:15])

for g in groups:
    print(g, top_genes[g])